import Utils.indexNode;
import com.fasterxml.jackson.core.JacksonException;
import com.fasterxml.jackson.databind.ObjectMapper;
import rtree.Rectangle;

import java.io.*;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class Silo {
    public indexNode dataLakeRoot;
    public ArrayList<indexNode> indexNodes = new ArrayList<indexNode>();
    public HashSet<indexNode> candidatesLeafnode = new HashSet<>();
    public static int a = 6;  //   3;   //
    public static int b = 7;  //   4;  //
    public static int c = 0;



    public void indexCheckinData(String path, int totalDatasetNumber) throws IOException,CloneNotSupportedException, ClassNotFoundException {
        long time = System.currentTimeMillis();
        generateSignatureFile(Main.datasetIdMapSignature, Main.datalakeRootNodes, indexNodes, Main.zcodeMBRCenter, Main.dimension,
        path, Main.minX, Main.minY, Main.rangeX, Main.rangeY, Main.resolution, totalDatasetNumber);
        long time11 = System.currentTimeMillis();
        indexNode parentNode = new indexNode(2);
        dataLakeRoot = indexDatasetKD(Main.datasetIdMapSignature, Main.datalakeRootNodes, parentNode, Main.dimension, Main.capacity,1);
        long time22 = System.currentTimeMillis();
        System.out.print("Create ball tree index = "+(time22- time11)+"ms, Number of datasets = "+ Main.datasetIdMapSignature.size()+". ");

    }

    public static indexNode indexDatasetKD(HashMap<String, ArrayList<Long>> datasetIdMapSignature, HashMap<String, indexNode> datalakeRootNodes, indexNode parentNode, int dimension, int leafSize, int upperGlobalId) {
        double minBox[] = new double[dimension];
        double maxBox[] = new double[dimension];
        double pivot[] = new double[dimension];
        for(int dim=0; dim<dimension; dim++) {
            minBox[dim] = Double.MAX_VALUE;
            maxBox[dim] = -1000000000;
        }
        // get the bounding box of all the nodes first,
        for(String i : datalakeRootNodes.keySet()) {
            indexNode aIndexNode = datalakeRootNodes.get(i);
            double nodeMax[] = aIndexNode.getMBRmax();
            double nodeMin[] = aIndexNode.getMBRmin();
            for(int dim=0; dim<dimension; dim++) {
                if(nodeMax[dim] > maxBox[dim])
                    maxBox[dim] = nodeMax[dim];//need to revise, as the radius is small due to weight
                if(nodeMin[dim] < minBox[dim])
                    minBox[dim] = nodeMin[dim];
                pivot[dim] += aIndexNode.getPivot()[dim];
            }
        }
        // get the pivot point, and the range in multiple dimension, and the radius
        int d=0;
        double maxrange = Double.MIN_VALUE;
        double radius = 0;
        for(int dim=0; dim<dimension; dim++) {
//            pivot[dim] = (maxBox[dim] + minBox[dim])/2;
            pivot[dim] = pivot[dim]/datalakeRootNodes.size();
//            radius += Math.pow((maxBox[dim] - minBox[dim])/2, 2);
            if((maxBox[dim] - minBox[dim]) > maxrange) {
                maxrange = maxBox[dim] - minBox[dim];
                d = dim;
            }
        }

        double leftUpBox[] = {minBox[0], maxBox[1]};
        double rightBottom[] = {maxBox[0], minBox[1]};

        double d1 = Math.max(distance2(pivot,minBox), distance2(pivot,leftUpBox));
        double d2 = Math.max(d1, distance2(pivot,rightBottom));
        double d3 = Math.max(d2, distance2(pivot, maxBox));
        radius = Math.sqrt(d3);
        //create a new leaf node and return
//        System.out.println("maxBox = "+maxBox[0]+", "+ maxBox[1]);
//        System.out.println("minBox = "+minBox[0]+", "+ minBox[1]);
//        System.out.println("pivot = "+pivot[0]+", "+pivot[1]);
//        System.out.println("radius =  "+radius);
//        System.out.println();
        indexNode a = new indexNode(dimension);
        a.setRadius(radius);// the radius need to be bigger.
        a.setPivot(pivot);
        a.setMBRmax(maxBox);
        a.setMBRmin(minBox);
        a.setUpperGlobalID(upperGlobalId);
        a.setParentNode(parentNode);
        if(datalakeRootNodes.size() <= leafSize) {
            HashMap<String, ArrayList<Long>> dataset = new HashMap<>();
            for(String i: datalakeRootNodes.keySet()) {
                indexNode root = datalakeRootNodes.get(i);
                ArrayList<Long> array = datasetIdMapSignature.get(i);
                dataset.put(i, array);
                a.addNodes(root);
//                System.out.println("id = "+root.getNodeid()+", radius = " + root.radius);
                a.addNodeIds(root.getNodeid());
                root.setParentNode(a);
            }
            //a.rootToDataset = -1
            a.setroot(-1);//label it as leaf node
            return a;
        }else {
            // dividing the space by the broadest dimension d
//                ArrayList<indexNode> rootleft = new ArrayList<indexNode>();
//                ArrayList<indexNode> rootright = new ArrayList<indexNode>();
            HashMap<String, indexNode> rootleft = new HashMap<>();
            HashMap<String, indexNode> rootright = new HashMap<>();
//            d=0;
            for(String i: datalakeRootNodes.keySet()) {
                indexNode aIndexNode = datalakeRootNodes.get(i);
                if(aIndexNode.getPivot()[d] < pivot[d]) {
                    rootleft.put(i, aIndexNode);
                }else {
                    rootright.put(i, aIndexNode);
                }
            }
            if(rootleft.isEmpty() || rootright.isEmpty()) {
                HashMap<String, ArrayList<Long>> dataset = new HashMap<>();
                for(String i: datalakeRootNodes.keySet()) {
                    indexNode root = datalakeRootNodes.get(i);
                    ArrayList<Long> array = datasetIdMapSignature.get(i);
                    dataset.put(i, array);
//                        System.out.println("id = "+root.getNodeid()+", radius = " + root.radius);
                    a.addNodes(root);
                    a.addNodeIds(root.getNodeid());
                    root.setParentNode(a);

                }
                a.setroot(-1);//label it as leaf node

            }else {
                a.addNodes(indexDatasetKD(datasetIdMapSignature, rootleft, a, dimension, leafSize, 2*upperGlobalId));
                a.addNodes(indexDatasetKD(datasetIdMapSignature, rootright, a, dimension, leafSize, 2*upperGlobalId+1));
                a.setroot(-2);//lable it as internal node
            }
            return a;
        }
    }

    public HashMap<Integer, HashSet<Integer>> geneOrloadGraphSer(String graphfilePath, double delta) throws JacksonException, IOException {
        HashMap<Integer, HashSet<Integer>> graphInfo33 = new HashMap<>();
//        ArrayList<ArrayList<Integer>> graphInfo33= new ArrayList<>();
        ObjectMapper objectMapper = new ObjectMapper();
        try {
            File file = new File(graphfilePath);
            if (file.exists()) {
//                System.out.println("graph is exists");
                graphInfo33 = deSerializationZcurve(graphfilePath);
//                FileInputStream inputStream = new FileInputStream(file);
//                int length = inputStream.available();
//                byte bytes[] = new byte[length];
//                inputStream.read(bytes);
//                inputStream.close();
//                String jsonString =new String(bytes, StandardCharsets.UTF_8);
//                graphInfo33 = objectMapper.readValue(jsonString, HashMap.class);
//                for (Integer s: graphString.keySet()){
//                    HashSet<Integer> hm = graphString.get(s);
////                    HashMap<String, Integer> hm = graphString.get(s);
//                    HashMap<String, Integer> hmInteger = new HashMap<String, Integer>();
//
//                    for (String s1: hm.keySet()){
//                        String key = s1;
//                        int value = (int)hm.get(s1);
//                        hmInteger.put(key , value);
//                    }
//                    graphInfo33.put(s, hmInteger);
//                }
            } else {
                for (String i : Main.datalakeRootNodes.keySet()) {
                    indexNode in = Main.datalakeRootNodes.get(i);
//                    ArrayList<Long> queryIn = datasetIdMapSignature.get(i);
//                    HashMap<String, Integer> edgeInfo = new HashMap<>();
                    HashSet<Integer> edgeInfo = new HashSet<>();
                    candidatesLeafnode.clear();
                    IBBranchAndBound(dataLakeRoot, in, delta); //     找候选集
                    for (indexNode candi : candidatesLeafnode) {
                        edgeInfo.add(Integer.parseInt(candi.getNodeid()));
//                        Set<indexNode> nodeList = candi.getNodelist();
//                        for (indexNode node: nodeList){
//                            if (distance(candi.getPivot(), in.getPivot())<= delta){
//                                if (!candi.getNodeid().equals(in.getNodeid())){
//                                    int frequence = 0;
//                                    edgeInfo.put(candi.getNodeid(), 0);
//                                }
//                            }
//                        }
                    }
                    graphInfo33.put(Integer.parseInt(i), edgeInfo);

                }
                System.out.println("success");
                SerializedZcurve(graphfilePath, graphInfo33);

//                System.out.println("isDiagonalSymmetric $$$$$$$$$$$$ = "+isDiagonalSymmetric(graphMatrix));
//                String jsonString1 = objectMapper.writerWithDefaultPrettyPrinter()
//                        .writeValueAsString(graphInfo33);
//                FileWriter fw = null;
//                try {
//                    fw = new FileWriter(graphfilePath);
//                    fw.write(jsonString1);
//                    fw.close();
////                    System.out.println("graph write successfully");
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
            }
        } catch(Exception e){
            e.printStackTrace();
        }
        return graphInfo33;
    }

    public static boolean isDiagonalSymmetric(int[][] matrix) {
        int n = matrix.length;
        boolean flag = true;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (matrix[i][j] != matrix[j][i]) {
                    flag = false;
                }
            }
        }
        return flag;
    }

    public HashMap<Integer, HashSet<Integer>> geneOrloadGraphJackson(String graphfilePath, double delta) throws JacksonException, IOException {
        HashMap<Integer, HashSet<Integer>> graphInfo33 = new HashMap<>();
//        ArrayList<ArrayList<Integer>> graphInfo33= new ArrayList<>();

        ObjectMapper objectMapper = new ObjectMapper();
        try {
            File file = new File(graphfilePath);
            if (file.exists()) {
//                System.out.println("graph is exists");
                FileInputStream inputStream = new FileInputStream(file);
                int length = inputStream.available();
                byte bytes[] = new byte[length];
                inputStream.read(bytes);
                inputStream.close();
                String jsonString =new String(bytes, StandardCharsets.UTF_8);
                HashMap<String, ArrayList<Integer>> graphString = objectMapper.readValue(jsonString, HashMap.class);
                for (String s: graphString.keySet()){
                    ArrayList<Integer> hm = graphString.get(s);
                    HashSet<Integer> hmInteger = new HashSet<Integer>();
                    for (Integer s1: hm){
                        hmInteger.add(s1);
                    }
                    graphInfo33.put(Integer.parseInt(s), hmInteger);
                }
            } else {


                for (String i : Main.datalakeRootNodes.keySet()) {
                    indexNode in = Main.datalakeRootNodes.get(i);
                    HashSet<Integer> edgeInfo = new HashSet<>();
                    candidatesLeafnode.clear();


                    IBBranchAndBound(dataLakeRoot, in, delta);

                    for (indexNode candi : candidatesLeafnode) {
                        edgeInfo.add(Integer.parseInt(candi.getNodeid()));
                    }

                    graphInfo33.put(Integer.parseInt(i), edgeInfo);
                }
                System.out.println("graph generate success");

                String jsonString1 = objectMapper.writerWithDefaultPrettyPrinter()
                        .writeValueAsString(graphInfo33);
                FileWriter fw = null;
                try {
                    fw = new FileWriter(graphfilePath);
                    fw.write(jsonString1);
                    fw.close();
                    System.out.println("graph write successfully");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } catch(Exception e){
            e.printStackTrace();
        }

        return graphInfo33;
    }


    public void getAllleafNode(indexNode root, HashSet<String> candidateList){
        if (root.getroot() == -1 ) {// leaf node
            Set<String> S =  root.getNodeIDList();
            for (String id: S){
//                if (datasetPrice.get(id)<= budget)
                    candidateList.add(id);
            }
        }else if(root.getroot() == -2){ //internal node
            Set<indexNode> listnode = root.getNodelist();
            for (indexNode aListNode: listnode){
                getAllleafNode(aListNode,candidateList);
            }
        }
    }
    public boolean pointContained(Rectangle rec, indexNode root, int dim){

        double[] mbrmin = rec.getLow().data;
        double[] mbrmax = rec.getHigh().data;
        ArrayList<Long> arrayList = Main.datasetIdMapSignature.get(root.getNodeid());
        for (long l: arrayList){
            double[] coordi = resolve(l, Main.resolution);
            for(int i=0; i<dim; i++){
                if (mbrmin[i] > coordi[i] || mbrmax[i] < coordi[i])
                    return false;
            }
        }
        return true;
    }

    public boolean contained(Rectangle rec, indexNode query, int dim){
        double[] mbrmin = rec.getLow().data;
        double[] mbrmax = rec.getHigh().data;
        double[] querymin = query.getMBRmin();
        double[] querymax = query.getMBRmax();
        for(int i=0; i<dim; i++){
            if (mbrmin[i] > querymin[i] || mbrmax[i] < querymax[i]){
                return false;
            }
        }
        return true;
//            if (mbrmin[0] <= querymin[0] && mbrmin[1] <= querymin[1] && mbrmax[0] >= querymax[0] && mbrmax[1] >= querymax[1]


    }
    public boolean intersected(Rectangle rec, indexNode query, int dim) {
        double[] mbrmin = rec.getLow().data;
        double[] mbrmax = rec.getHigh().data;
        double[] querymin = query.getMBRmin();
        double[] querymax = query.getMBRmax();
        for(int i=0; i<dim; i++){
            if(mbrmin[i] > querymax[i] || mbrmax[i] < querymin[i]) {
                return false;
            }
        }

        return true;
    }


    public void IBBranchAndBound(indexNode root, indexNode query, double delta){
        if (root.getroot() == -1 ) { // leaf node
            Set<indexNode> nodeList = root.getNodelist();
            for (indexNode node: nodeList) {
                if (distance(node.getPivot(), query.getPivot()) <= delta) {
                    if (node.getNodeid() != query.getNodeid()) {
                        candidatesLeafnode.add(node);
                    }
                }
            }
        }else if(root.getroot() == -2){ //internal node
            Set<indexNode> listnode = root.getNodelist();
            for (indexNode aListNode: listnode){
                if (distance(aListNode.getPivot(), query.getPivot()) - aListNode.getRadius() >delta){
                    //filter
                }else { //if (distance(aListNode.getPivot(), query.getPivot()) - aListNode.getRadius() <= delta)
                    IBBranchAndBound(aListNode,query, delta);
                }
            }
        }

    }

    public void generateSignatureFile(
            HashMap<String, ArrayList<Long>> datasetIdMapSignature, HashMap<String, indexNode> datalakeRootNodes,
            ArrayList<indexNode> indexNodes, HashMap<String, double[][]> zcodeMBRCenter, int dimension, String Path,
            double minX, double minY, double rangeX, double rangeY, int resolution, int totalDatasetNumber) throws IOException, CloneNotSupportedException {
        String signaturePath = Main.preprocessFilePath+ Main.siloName +"/";
        File f =new File(signaturePath);
        if  (!f.exists()  && !f.isDirectory()) {
            f.mkdir();
            System.out.println("create success");
        }
//        HashMap<String, HashMap<Long, Double>> zcodeMap = new HashMap<String, HashMap<Long, Double>>();
        String filePath = signaturePath  + "resolution-"+resolution + "-number"+totalDatasetNumber + ".txt";
//        resolution +"-delta"+delta+"-number"+totalDatasetNumber
        try {
            File file = new File(filePath);
            long time1 = System.currentTimeMillis();
            if (file.exists()) {
                //直接读文件
//                System.out.println("file = "+filePath+" is exist!!!");
                try(BufferedReader br = new BufferedReader(new FileReader(filePath))){
                    String strLine;
                    while ((strLine = br.readLine()) != null) {
                        HashMap<Long, Double> zcodeHashMap = new HashMap<>();
                        ArrayList<Long> arrayList = new ArrayList<>();
                        String[] splitString = strLine.split(";");
//                        int id = Integer.parseInt(splitString[0]);
                        String[] idAndMbr = splitString[0].split(",");
                        String mappedStringId = idAndMbr[0];
                        String datasetId = idAndMbr[1];
                        double[][] MBRs = new double[3][2];

                        MBRs[0][0] = Double.parseDouble(idAndMbr[2]);
                        MBRs[0][1] = Double.parseDouble(idAndMbr[3]);
                        MBRs[1][0] = Double.parseDouble(idAndMbr[4]);
                        MBRs[1][1] = Double.parseDouble(idAndMbr[5]);
                        MBRs[2][0] = Double.parseDouble(idAndMbr[6]);
                        MBRs[2][1] = Double.parseDouble(idAndMbr[7]);

                        Main.integerMaptoDatasetId.put(Integer.parseInt(idAndMbr[0]), datasetId);

                        double max[], min[];
                        max = MBRs[0];
                        min = MBRs[1];

                        for (int i = 1; i<splitString.length; i++){
                            arrayList.add(Long.parseLong(splitString[i]));
                        }

                        double pivot[] = MBRs[2];
                        double radius = Math.max(distance(pivot, MBRs[0]), distance(pivot, MBRs[1]) );

                        radius = Math.sqrt(radius);
                        indexNode rootNode = new indexNode(dimension);
                        rootNode.setMBRmax(max);
                        rootNode.setMBRmin(min);
                        rootNode.setRadius(radius);
                        rootNode.setPivot(pivot);
//                        System.out.println("setNodeid3 = "+ id);
                        rootNode.setNodeid(mappedStringId);
//                        rootNode.setDeviceID(siloID);
//                        zcodeMap.put(numericId, zcodeHashMap);
                        zcodeMBRCenter.put(mappedStringId, MBRs);
                        datasetIdMapSignature.put(mappedStringId, arrayList);
                        datalakeRootNodes.put(mappedStringId,rootNode);
                        indexNodes.add(rootNode);
                    }
                }catch (IOException e) {
                    e.printStackTrace();
                }
            }else {
                File[] listFiles = new File(Path).listFiles();
                if (totalDatasetNumber<listFiles.length){
                    for (int t=0; t<totalDatasetNumber; t++) {
                        File filename = listFiles[t];
                        String datasetid = filename.getName();
                        String dataPath = Path + datasetid;
                        int integerId = t;
                        Main.integerMaptoDatasetId.put(integerId, datasetid);
                        double[][] MBRCenter = new double[3][2];
                        HashSet<Long> hs = new HashSet<>();
                        HashMap<Long, Double> zcodeHashMap = generateSignature(hs, MBRCenter, minX, minY, rangeX, rangeY, resolution, dataPath);
                        if (zcodeHashMap.size() > 0) {
                            ArrayList<Long> arrayList = new ArrayList<>(hs);
                            double max[], min[], pivot[], radius = 0;
                            max = MBRCenter[0];
                            min = MBRCenter[1];
                            pivot = MBRCenter[2];
                            radius = Math.max(distance(pivot, max), distance(pivot, min));
                            indexNode rootNode = new indexNode(dimension);
                            rootNode.setMBRmax(max);
                            rootNode.setMBRmin(min);
                            rootNode.setRadius(radius);
                            rootNode.setPivot(pivot);
                            String mappedStringId = Integer.toString(integerId);
                            rootNode.setNodeid(mappedStringId);
//                            zcodeMap.put(stringId, zcodeHashMap); //这里zcodeHashMap存储了cell id 和 number of points
                            zcodeMBRCenter.put(mappedStringId, MBRCenter);
                            datasetIdMapSignature.put(mappedStringId, arrayList);
                            datalakeRootNodes.put(mappedStringId, rootNode);
                            indexNodes.add(rootNode);
                        }
                    }
                }else {
                    totalDatasetNumber = listFiles.length;
//                    File filename: listFiles
                    for (int t = 0; t < totalDatasetNumber; t++) {
                        File filename = listFiles[t];
                        String datasetid = filename.getName();
                        String dataPath = Path + datasetid;
                        int integerId = t;
                        Main.integerMaptoDatasetId.put(integerId, datasetid);
                        double[][] MBRCenter = new double[3][2];
                        HashSet<Long> hs = new HashSet<>();
                        HashMap<Long, Double> zcodeHashMap = generateSignature(hs, MBRCenter, minX, minY, rangeX, rangeY, resolution, dataPath);
                        HashMap<Integer, HashSet<Integer>> sliceInformation = new HashMap<>();
                        if (zcodeHashMap.size() > 0) {
                            ArrayList<Long> arrayList = new ArrayList<>(hs);
                            double max[], min[], pivot[];
                            max = MBRCenter[0];
                            min = MBRCenter[1];
                            pivot = MBRCenter[2];

                            double radius = Math.max(distance(pivot, max), distance(pivot, min));
                            indexNode rootNode = new indexNode(dimension);
                            rootNode.setMBRmax(max);
                            rootNode.setMBRmin(min);
                            rootNode.setRadius(radius);
                            rootNode.setPivot(pivot);
//                            System.out.println("setNodeid1 = "+ datasetid);
                            String mappedStringId = Integer.toString(integerId);
                            rootNode.setNodeid(mappedStringId);
//                            zcodeMap.put(numericId, zcodeHashMap); //这里zcodeHashMap存储了cell id 和 number of points
                            zcodeMBRCenter.put(mappedStringId, MBRCenter);
                            datasetIdMapSignature.put(mappedStringId, arrayList);
                            datalakeRootNodes.put(mappedStringId, rootNode);
                            indexNodes.add(rootNode);
                        }
                    }
                }

                writeFileAndMBR(zcodeMBRCenter, Main.integerMaptoDatasetId, filePath, datasetIdMapSignature);
            }
            long time2 = System.currentTimeMillis();

        }catch(Exception e){
            e.printStackTrace();
        }
//        return zcodeMap;
    }


    public static double[] resolve(Long key, int resolution){
        int d0 = (int)(key % Math.pow(2,resolution)); //x
        int d1 = (int)(key / Math.pow(2,resolution));  //y
        double[] d = {d0, d1};
        return d;
    }
    public static long combine(double[] coordi, int resolution){
        long d = 0;
        double t = Math.pow(2,resolution);
        d = (long)(coordi[1]*t+coordi[0]);
        return d;
    }


    public static void writeFile(String filePath, HashMap<Integer, HashMap<Long,Double>> dataSetMap) throws IOException {
        FileWriter fw = null;
        try {
            File file = new File(filePath);
            file.createNewFile();
            fw = new FileWriter(filePath);
            for (int id : dataSetMap.keySet()) {
                HashMap<Long,Double> map1 = dataSetMap.get(id);
                fw.write(String.valueOf(id));
                fw.write(";");
                for (long key: map1.keySet()){
                    fw.write(String.valueOf(key));
                    fw.write(",");
                    String aa = String.valueOf(map1.get(key));
                    BigDecimal db = new BigDecimal(aa);
                    String ii = db.toPlainString();
                    fw.write(ii+";");
                }
                fw.write("\n");
//                    System.out.println();
            }
            fw.close();

        }catch(Exception e){
            e.printStackTrace();
        }
    }

    public  HashMap<Long,Double> generateSignature(HashSet<Long> hs, double[][] MBRCenter, double minX, double minY, double rangeX, double rangeY, int resolution, String Path) throws CloneNotSupportedException{
        double blockSizeX = rangeX/Math.pow(2,resolution);
        double blockSizeY = rangeY/Math.pow(2,resolution);
        HashMap<Long,Integer> ZorderHis = new HashMap<>();
        HashMap<Long,Double> ZorderDensityHist = new HashMap<>();
        long t = (long)Math.pow(2,resolution);
        if (Main.siloName.startsWith("argoverse")){
            a = 6;
            b = 7;
            c = 2;
        }else if (Main.siloName.startsWith("baidu")){
            a = 1; //lng
            b = 2; //lat
        }else{
            a = 0;
            b = 1;
        }
        double[] maxMBR = {-1000000000, -1000000000};
        double[] minMBR = {Double.MAX_VALUE, Double.MAX_VALUE};
        try {
            String record = "";
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(Path)));
            reader.readLine();
            while((record=reader.readLine())!=null) {
                String[] fields = record.split(",");
                double x_coordinate = Double.parseDouble(fields[a]);
                double y_coordinate = Double.parseDouble(fields[b]);
                if (fields[a].matches("-?\\d+(\\.\\d+)?")){ // && x_coordinate>-180 && x_coordinate < 180 && y_coordinate> -90 && y_coordinate <90

                    double x = x_coordinate - minX;
                    double y = y_coordinate - minY;
                    long X = (long) (x / blockSizeX);  //行row
                    long Y = (long) (y / blockSizeY);  //列 col
                    long id = Y * t + X;

                    hs.add(id);
                    Integer num = ZorderHis.getOrDefault(id, new Integer(0));
                    ZorderHis.put(id, num + 1);
                    ZorderDensityHist.put(id, (double) (num + 1));
                    if (maxMBR[0] < X)
                        maxMBR[0] = X;
                    if (maxMBR[1] < Y)
                        maxMBR[1] = Y;
                    if (minMBR[0] > X)
                        minMBR[0] = X;
                    if (minMBR[1] > Y)
                        minMBR[1] = Y;
                }
            }


            double[] temp = new double[2];
            for (long l:hs){
                double[] coor = resolve(l, resolution);
                temp[0] += coor[0];
                temp[1] += coor[1];
            }

            double[] center = {temp[0]/hs.size(), temp[1]/hs.size()};
            MBRCenter[0] = maxMBR.clone();
            MBRCenter[1] = minMBR.clone();
            MBRCenter[2] = center;

//            System.out.println(" MBR(gen) = "+maxMBR[0]+", "+maxMBR[1]+", "+minMBR[0]+", "+minMBR[1] +", "+center[0]+", "+center[1] );
            reader.close();
        }catch (IOException e) {
//            e.printStackTrace();
//            System.out.println("this dataset name isn't exist");
        }

        return ZorderDensityHist;
    }


    public static void SerializedZcurve(String file, HashMap<Integer, HashSet<Integer>> result) {
        try {
            FileOutputStream fos = new FileOutputStream(file);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(result);
            oos.close();
            fos.close();
            System.out.println("Serialized result HashMap data is saved in hashmap.ser");
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public static HashMap<Integer, HashSet<Integer>> deSerializationZcurve(String file) {

        HashMap<Integer, HashSet<Integer>> result;
        try {
            FileInputStream fis = new FileInputStream(file);
            ObjectInputStream ois = new ObjectInputStream(fis);
            result = (HashMap) ois.readObject();
            ois.close();
            fis.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
            return null;
        } catch (ClassNotFoundException c) {
            System.out.println("Class not found");
            c.printStackTrace();
            return null;
        }
        return result;
    }

    public static List<indexNode> getAllIndexNode(File file) throws FileNotFoundException, IOException, ClassNotFoundException {
        FileInputStream fis = new FileInputStream(file);
        ObjectInputStream ois = new ObjectInputStream(fis);
        List<indexNode> objArr = new ArrayList<indexNode>();
        indexNode p = null;
        while (fis.available() > 0) {
            p = (indexNode) ois.readObject();
            objArr.add(p);
        }
        ois.close();
        return objArr;
    }
    public static String bytesToString(byte[] bytes) {
        //转换成base64
        return org.apache.commons.codec.binary.Base64.encodeBase64String(bytes);
    }

    public static double distance2(double[] x, double[] y) {
        double d = 0.0;
        for (int i = 0; i < x.length; i++) {
            d += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return d;
    }
    public static double distance(double[] x, double[] y) {
        return Math.sqrt(distance2(x, y));
    }

    public static void writeFileAndMBR(HashMap<String, double[][]> zcodeMBRCenter, HashMap<Integer, String> integerMaptoDatasetId, String filePath, HashMap<String, ArrayList<Long>> dataSetMap) throws IOException {
//        System.out.println("dataSetMap.size() = "+dataSetMap.size());

        FileWriter fw = null;
        try {
            File file = new File(filePath);
            file.createNewFile();
            fw = new FileWriter(filePath);
            int count = 0;
            for (String id : dataSetMap.keySet()) {
                count++;
//                System.out.println(id+", "+integerMaptoDatasetId.get(Integer.parseInt(id)));
                ArrayList<Long> map1 = dataSetMap.get(id);
                fw.write(id+",");
                fw.write(integerMaptoDatasetId.get(Integer.parseInt(id))+",");
                double[][] MBRs = zcodeMBRCenter.get(id);
                fw.write(MBRs[0][0]+","); // leftBottom
                fw.write(MBRs[0][1]+",");
                fw.write(MBRs[1][0]+","); // rightUp
                fw.write(MBRs[1][1]+",");
                fw.write(MBRs[2][0]+","); // center coordinate
                fw.write(MBRs[2][1]+";");
//                fw.write(";");
                for (long key: map1){
                    fw.write(String.valueOf(key));
                    fw.write(";");
//                    fw.write(",");
//                    String aa = String.valueOf(map1.get(key));
//                    BigDecimal db = new BigDecimal(aa);
//                    String ii = db.toPlainString();
//                    fw.write(ii+";");
                }
                fw.write("\n");

            }

            fw.close();
            System.out.println("count = "+count +", resolution.txt write successfully");

        }catch(Exception e){
            e.printStackTrace();
        }
    }


}

