import Utils.indexNode;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
public class Main {

    public static double minX = -180;
    public static double minY = -90;
    public static double rangeX = 360;
    public static double rangeY = 180;
    public static int dimension = 2;
    public static int capacity = 10;
    public static HashMap<String, ArrayList<Long>> datasetIdMapSignature = new HashMap<>();
    public static HashMap<String, double[][]> zcodeMBRCenter = new HashMap<>();

    public static String directory= "D:\\GitHub\\SpatialDataAcquisition\\src\\Data\\DataCollections\\";
    public static String preprocessFilePath = "D:\\GitHub\\SpatialDataAcquisition\\src\\Data\\";
    public static String siloName = "openGeo";
    public static String outFilePath1 = "D:\\GitHub\\SpatialDataAcquisition\\src\\Data\\result\\";
    public static int totalDatasetNumber = 1000;
    public static int resolution = 11;
    public static double delta = 0;
    public static double budget = 0.1;
    public static String graphfilePath = "";
    public static Silo silo;
    public static double p_min = 1000000000000000000.0;
    public static double p_max = 0;
    public static indexNode BFS_root = new indexNode(2);
    public static HashMap<Integer, HashSet<Integer>> graphInfo = new HashMap<>();
    public static HashMap<String, HashSet<String>> queryWindowGraphInfo = new HashMap<>();
    public static HashMap<String, indexNode> datalakeRootNodes = new HashMap<>();
    public static HashMap<String, Double> datasetPrice = new HashMap<>();

    public static HashMap<Integer, String> integerMaptoDatasetId = new HashMap<>();

    public static void main(String[] args) throws IOException, CloneNotSupportedException, ClassNotFoundException{

        boolean wirteLog = false;
        String[] siloNameList = {"trackable"}; // trackable, public, identifiable, btaa-clean, openGeo
        int[] resolutionList = {11};
        double[] deltaList = {10};
        double[] budgetList = {0.001};
        String[] algorithmList = {"DSA", "CMC+MC", "CMC+MG", "DPSA+BA", "BGPA"};
        for (int t=0; t <siloNameList.length; t++){
            siloName = siloNameList[t];
            for (int i = 0; i<resolutionList.length; i++) {
                datasetIdMapSignature = new HashMap<>();
                datalakeRootNodes = new HashMap<>();
                resolution = resolutionList[i];
                datasetPrice = new HashMap<>();
                String index = "IBtree"; // invertedIndex // IBtree
                preProcessing(index);

                double totalPrice = assignPrice();

                for (int j = 0; j < deltaList.length; j++) {
                    delta = deltaList[j];
                    graphfilePath = preprocessFilePath + siloName + "/" + "graph-" + "reso" + resolution + "-delta" + delta + "-number" + totalDatasetNumber + ".txt"; //ser
                    graphInfo = new HashMap<>();

                    //load graph
                    long time11 = System.currentTimeMillis();
                    graphInfo = silo.geneOrloadGraphJackson(graphfilePath, delta);
                    long time12 = System.currentTimeMillis();

                    System.out.println("graphInfo.keySet().size() = "+graphInfo.keySet().size());
                    for (int s1 : graphInfo.keySet()) {
                        HashSet<Integer> edgeInfo = graphInfo.get(s1);
//                            System.out.print(s1+", price="+datasetPrice.get(s1)+", "+edgeInfo.size());
                        HashSet<String> subEdgeInfo = new HashSet<>();
                        for (Integer s2 : edgeInfo) {
                            subEdgeInfo.add(String.valueOf(s2));
                        }
                        queryWindowGraphInfo.put(String.valueOf(s1), subEdgeInfo);
                    }

                    HashMap<String, HashSet<String>> subGraphNodeIdLists = new HashMap<>();
                    HashSet<String> allNodes = new HashSet<>();
                    for (String s1 : queryWindowGraphInfo.keySet()) {
                        allNodes.add(s1);
                    }
                    for (String s1 : queryWindowGraphInfo.keySet()) {
                        if (allNodes.contains(s1)) {
                            HashSet<String> hs = new HashSet<>();
                            hs.add(s1);

                            Queue<String> queue = new LinkedList<>();
                            queue.add(s1);
                            allNodes.remove(s1);
                            while (!queue.isEmpty()) {
                                String id = queue.poll();
                                HashSet<String> edgeInfo = queryWindowGraphInfo.get(id);
                                for (String s2 : edgeInfo) {
                                    if (!hs.contains(s2)) {
                                        hs.add(s2);
                                        queue.add(s2);
                                        allNodes.remove(s2);
                                    }
                                }
                            }
                            subGraphNodeIdLists.put(s1, hs);
                        }
                    }

                    PriorityQueue<relaxIndexNode> subGraphNodeIdListCovandPrice = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
                        @Override
                        public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                            return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
                        }
                    }); // large -> small
                    for (String s1 : subGraphNodeIdLists.keySet()) {
                        double coverageNumber = 0;
                        double price = 0;
                        HashSet<String> hs = subGraphNodeIdLists.get(s1);
                        for (String s2: hs){
                            coverageNumber += datasetIdMapSignature.get(s2).size();
                            price += datasetPrice.get(s2);
                        }
                        subGraphNodeIdListCovandPrice.add(new relaxIndexNode(s1, coverageNumber, price));
                    }

                    HashMap<String, HashMap<String, HashSet<String>>> subGraphS = new HashMap<>();
                    HashMap<String, HashMap<String, indexNode>> subGraphIndexNodeS = new HashMap<>();
                    for (String s : subGraphNodeIdLists.keySet()) {
                        HashSet<String> subGraphNodeIds = subGraphNodeIdLists.get(s);
                        HashMap<String, HashSet<String>> subGraph = new HashMap<>();
                        HashMap<String, indexNode> subGraphIndexNode = new HashMap<>();
                        for (String s1 : subGraphNodeIds) {
                            HashSet<String> edgeInfo = queryWindowGraphInfo.get(s1);
                            HashSet<String> subEdgeInfo = new HashSet<>();
                            for (String s2 : edgeInfo) {
                                if (subGraphNodeIds.contains(s2)) {
                                    subEdgeInfo.add(s2);
                                }
                            }
                            subGraph.put(s1, subEdgeInfo);
                            subGraphIndexNode.put(s1, new indexNode(2));
                        }
                        subGraphS.put(s, subGraph);
                        subGraphIndexNodeS.put(s, subGraphIndexNode);
                    }


//                        statisticInfo(graphInfo, subGraphS);



                    for (int k = 0; k < budgetList.length; k++) {
                        budget = budgetList[k];   //  + "Number"+totalDatasetNumber
                        String outFilePath2 = outFilePath1+ siloName+"/"  + "1227-reso" + resolution +"-delta"+delta+"-budget"+budget+"-number"+totalDatasetNumber+".txt";
                        budget = budgetList[k] * totalPrice; // 100;  //180 ;//budgetList[k] * totalPrice;
                        File file = new File(outFilePath1 + siloName+"/");
                        if (!file.exists()){
                            file.mkdir();
                            System.out.println("create success 1");
                        }

                        File file2 = new File(outFilePath2);

                        if (!file2.exists()) {
                            if (wirteLog) {
                                file2.createNewFile();
                                System.out.println("create success 2");
                                PrintStream out = new PrintStream(outFilePath2);
                                System.setOut(out);
                            }
                            System.out.print("Load graph time = "+ (time12 - time11)+" ms, ");
                            System.out.println("subGraph.size() = " + subGraphNodeIdLists.size());
                            System.out.println("siloName = "+siloName+",  "+ "resolution = "+resolution+",  "+ "Budget = "+budget/totalPrice+ ", Delta = "+delta);
                            System.out.println(siloName+"'s total Price = "+totalPrice);
                            System.out.println("datasetIdMapSignature.size() = "+datasetIdMapSignature.size()+",  queryWindowGraphInfo.size() = " + queryWindowGraphInfo.size() + ",  subGraph.size() = " + subGraphNodeIdLists.size());


                            for (int al= 0; al<algorithmList.length; al++){
                                String algo = algorithmList[al];

                                long time31 = 0;
                                long time32 = 0;
                                long time = 0;
                                double totalCoverage = 0.0;
                                int iteration3 = 0;
                                int prune = 0;
                                HashMap<String, LinkedHashMap<String, Integer>> IBGPAResultS = new HashMap<>();
                                ArrayList<Long> findCenterTime = new ArrayList<>();
                                int IBGPAMaxCoverage;
                                int coverageNumber;
                                HashMap<String, Integer> IBGPAResultSHM = new HashMap<>();

                                switch (algo){
                                    case "visual-random":
                                        System.out.println("visual-random");
                                        List<String> keyList = new ArrayList<>(datasetIdMapSignature.keySet());
                                        HashSet<String> resultPath = new HashSet<>();
                                        boolean flag = true;
                                        int count = 0;
                                        while (getCurrentCost(resultPath) < budget && flag == true){
//                                            Random random = new Random(); // 创建随机数生成器
//                                            int randomIndex = random.nextInt(keyList.size()); // 随机选择一个键
                                            int randomIndex = count;
                                            String randomKey = keyList.get(randomIndex);
                                            if (!resultPath.contains(randomKey)){
                                                if (getCurrentCost(resultPath)+datasetPrice.get(randomKey)<= budget){
                                                    resultPath.add(randomKey);
                                                }else {
                                                    flag = false;
                                                }
                                            }
                                            count++;
                                        }
                                        HashSet<Long> coveredCells = new HashSet<>();

                                        for (String s : resultPath) {
                                            int tt = coveredCells.size();
                                            coveredCells.addAll(datasetIdMapSignature.get(s));
//                                            System.out.println(s+" ,"+integerMaptoDatasetId.get(Integer.parseInt(s))+", coverage incremental =" +  (coveredCells.size()-tt)+"; ");
                                        }
                                        System.out.println("random selected result, C.size = "+coveredCells.size());
                                        break;
                                    case "DSA":
                                        System.out.println("###############################");
                                        System.out.println(" begin DSA" + ",  budget = " + budget);
                                        long time1 = System.currentTimeMillis();
                                        LinkedHashMap<String, Integer> result = simpleGreedy2(graphInfo);
                                        long time2 = System.currentTimeMillis();
                                        System.out.println(", runningTime = " + (time2 - time1) + "ms");
                                        System.out.println();
                                        break;

                                    case "DPSA+BA":
                                        System.out.println("###############################");
                                        System.out.println(" begin DPSA+BA," + " budget = " + budget);
                                        time31 = System.currentTimeMillis();
                                        //                                    HashMap<String, LinkedHashMap<String, Integer>> IBGPAResultS = new HashMap<>();
                                        //                                    ArrayList<Long> findCenterTime = new ArrayList<>();
                                        iteration3 = 0;
                                        for (String s : subGraphNodeIdLists.keySet()) { //relaxIndexNode re: getDeletedSubgraphNodes(subGraphNodeIdListCovandPrice)
//                                            String s = re.resultId;
                                            getPminPmax(subGraphNodeIdLists.get(s));
                                            iteration3++;
                                            HashMap<String, HashSet<String>> subGraph = subGraphS.get(s);
                                            HashMap<String, indexNode> subGraphIndexNode = subGraphIndexNodeS.get(s);
                                            LinkedHashMap<String, Integer> rootLeafIncreHm = new LinkedHashMap<>();
                                            if (subGraph.size() > 1) {
//                                                System.out.println("subgrah-" + iteration3 +", "+ integerMaptoDatasetId.get(Integer.parseInt(s)) + ": pmax= " + outputPrice(p_max) + ", pmin= " + outputPrice(p_min) + ", ");
                                                rootLeafIncreHm = IBGPA(subGraphIndexNode, subGraph, true, findCenterTime);
                                            } else {
                                                for (String center : subGraph.keySet()) {
                                                    if (datasetPrice.get(center)<=budget)
                                                        rootLeafIncreHm.put(center, datasetIdMapSignature.get(center).size());
//                                                        System.out.println("subgrah(single)- "+iteration3 +" = "+center +", root node coverage = "+ datasetIdMapSignature.get(center).size() +", getCurrentCost(resultPath) = "+datasetPrice.get(center));
                                                }
                                            }
                                            IBGPAResultS.put(s, rootLeafIncreHm);
                                        }

                                        IBGPAMaxCoverage = findBest(IBGPAResultS);
                                        time32 = System.currentTimeMillis();

                                        time = 0;
                                        for (long l:findCenterTime){
                                            time += l;
                                        }
                                        System.out.println("BGPA+ BFS-based acceleration:" + IBGPAMaxCoverage  + ", find center time = " + time + "ms" + ", running time = " + (time32 - time31) + "ms");
                                        System.out.println();
                                        break;


                                    case "BGPA":
                                        System.out.println("###############################");
                                        System.out.println(" begin BGPA," + " budget = " + budget);
                                        time31 = System.currentTimeMillis();
                                        //
                                        iteration3 = 0;
                                        for (String s : subGraphNodeIdLists.keySet()) { //
                                            //            System.out.println();
                                            getPminPmax(subGraphNodeIdLists.get(s));
                                            iteration3++;
                                            HashMap<String, HashSet<String>> subGraph = subGraphS.get(s);
                                            HashMap<String, indexNode> subGraphIndexNode = subGraphIndexNodeS.get(s);
                                            LinkedHashMap<String, Integer> rootLeafIncreHm = new LinkedHashMap<>();
                                            if (subGraph.size() > 1) {
//                                            System.out.print("subgrah-" + iteration3 +", "+ integerMaptoDatasetId.get(Integer.parseInt(s)) +", subGraph.size() = " + subGraph.size() +": pmax= " + outputPrice(p_max) + ", pmin= " + outputPrice(p_min) + ", ");
                                                rootLeafIncreHm = IBGPA(subGraphIndexNode, subGraph, false, findCenterTime);
//                                            System.out.println();
                                            } else {
                                                //                System.out.println("this graph only has one node");
                                                for (String center : subGraph.keySet()) {
                                                    if (datasetPrice.get(center)<=budget){
                                                        rootLeafIncreHm.put(center, datasetIdMapSignature.get(center).size());
//                                                    System.out.println("iteration! "+iteration3 +" = "+center +", root node coverage = "+ datasetIdMapSignature.get(center).size() +", getCurrentCost(resultPath) = "+datasetPrice.get(center));
                                                    }
                                                }
                                            }
                                            IBGPAResultS.put(s, rootLeafIncreHm);
                                            //                                    break;
                                        }
                                        IBGPAMaxCoverage = findBest(IBGPAResultS);
                                        time32 = System.currentTimeMillis();
                                        time = 0;
                                        for (long l:findCenterTime){
                                            time += l;
                                        }
                                        System.out.println("IBGPA:" + IBGPAMaxCoverage + ", find center time = " + time + "ms" + ", running time = " + (time32 - time31) + "ms");
                                        System.out.println();
                                        break;
                                    case "CMC+MG":
                                        // for each subgraph find the result;
                                        System.out.println("###############################");
                                        System.out.println(" begin CMC + MG" + ",  budget = " + budget);
                                        time31 = System.currentTimeMillis();
                                        totalCoverage = 0.0;
                                        iteration3 = 0;
                                        prune = 0;
                                        for (String s: subGraphNodeIdLists.keySet()) {
                                            iteration3++;
                                            HashMap<String, HashSet<String>> subGraph = subGraphS.get(s);
                                            HashMap<String, indexNode> subGraphIndexNode = subGraphIndexNodeS.get(s);
                                            coverageNumber = 0;
                                            getPminPmax(subGraphNodeIdLists.get(s));

                                            if (subGraph.size() > 1) {
                                                //                                        System.out.print("subgrah-" + iteration3 + ": pmax= " + outputPrice(p_max) + ", pmin= " + outputPrice(p_min) + ", ");
                                                coverageNumber = CombinatorialAlgorithmSpeedMaxGain(subGraphIndexNode, subGraph);
                                                //                                            System.out.println("subgrah-" + iteration3 +": coverage= " +coverageNumber);
                                            } else {
                                                //                System.out.println("this graph only has one node");
                                                for (String center : subGraph.keySet()) {
                                                    if (datasetPrice.get(center)<=budget)
                                                        coverageNumber =  datasetIdMapSignature.get(center).size();
                                                    //                                                System.out.println("subgrah!-" + iteration3 + ": coverage= " +coverageNumber);
                                                    //                    System.out.println("iteration! "+1 +" = "+center +", root node coverage = "+ datasetIdMapSignature.get(center).size() +", getCurrentCost(resultPath) = "+datasetPrice.get(center));
                                                }
                                            }
                                            IBGPAResultSHM.put(s, coverageNumber);
                                            if (coverageNumber> totalCoverage)
                                                totalCoverage = coverageNumber;
                                        }
                                        //                                int IBGPAMaxCoverage =  findBest(IBGPAResultS);
                                        time32 = System.currentTimeMillis();
                                        System.out.println("CMC + MG:" + totalCoverage + ", running time = " + (time32 - time31) + "ms");
                                        System.out.println();
                                        break;

                                    case "CMC+MC":
                                        // //                           for each subgraph find the result;
                                        System.out.println("###############################");
                                        System.out.println(" begin CMC + MC" + ",  budget = " + budget);
                                        time31 = System.currentTimeMillis();
                                        totalCoverage = 0.0;
                                        iteration3 = 0;
                                        totalCoverage = 0;
//                                        System.out.println("subGraphNodeIdLists.keySet().size() = "+subGraphNodeIdLists.keySet().size());
                                        for (String s: subGraphNodeIdLists.keySet()) {
                                            iteration3++;
                                            HashMap<String, HashSet<String>> subGraph = subGraphS.get(s);
                                            HashMap<String, indexNode> subGraphIndexNode = subGraphIndexNodeS.get(s);
                                            coverageNumber = 0;
                                            getPminPmax(subGraphNodeIdLists.get(s));
                                            if (subGraph.size() > 1) {
                                                coverageNumber =  CombinatorialAlgorithmSpeedMaxSize(subGraphIndexNode, subGraph);
                                            }else {
                                                //                System.out.println("this graph only has one node");
                                                for (String center : subGraph.keySet()) {
                                                    if (datasetPrice.get(center)<=budget)
                                                        coverageNumber =  datasetIdMapSignature.get(center).size();
                                                }
                                            }
//                                            System.out.println("coverageNumber = "+ coverageNumber );
                                            IBGPAResultSHM.put(s, coverageNumber);
                                            if (coverageNumber> totalCoverage)
                                                totalCoverage = coverageNumber;
                                        }

                                        time32 = System.currentTimeMillis();
                                        System.out.println("CMC + MC: " + totalCoverage + ", running time = " + (time32 - time31) + "ms");
                                        System.out.println();
                                        break;

                                    default:
                                        System.out.println("algorithm name is error");
                                        break;

                                }
                            }
                        }
                    }

//*/

                }
            }
        }


        System.out.println("End!!!!");
    }

    public static void statisticInfo(HashMap<Integer, HashSet<Integer>> graphInfo, HashMap<String, HashMap<String, HashSet<String>>> subGraphS){
        // 1,统计边的数量
        int edgeNumber = 0;
        for (int i: graphInfo.keySet()){
            HashSet<Integer> edgeInfo = graphInfo.get(i);
            edgeNumber += edgeInfo.size();
        }
        System.out.println("edge Number = "+edgeNumber/2);

        //  2,average degree
        double AvgDegree = (double) edgeNumber/graphInfo.size();
        System.out.println("Average degree = "+AvgDegree);


        //  3,number of subgraphs
        System.out.println("number of subgraph = "+subGraphS.size());

        //  4，number of average cells
        int numberOfCells = 0;
        for (String s: datasetIdMapSignature.keySet()){
            numberOfCells+= datasetIdMapSignature.get(s).size();
        }
        double avgCells = numberOfCells/datasetIdMapSignature.size() ;
        System.out.println("number of average cells  = "+ avgCells);
    }




    public static void preProcessing(String index) throws IOException,CloneNotSupportedException, ClassNotFoundException{
        silo = new Silo();
        silo.indexCheckinData(directory+ siloName+"/", totalDatasetNumber);
    }

    public static LinkedHashMap<String, Integer> simpleGreedy2(HashMap<Integer, HashSet<Integer>> subGraph){

        LinkedHashMap<String, Integer> resultHs1 = new LinkedHashMap<>();
        int H1 = 0;
        H1 = simpleGreedySubFunction2(resultHs1, subGraph, false);
        int H2 = 0;
        LinkedHashMap<String, Integer> resultHs2 = new LinkedHashMap<>();
        resultHs2 = new LinkedHashMap<>();
//        H2 = simpleGreedySubFunction2(resultHs2, subGraph, true);
        if (H1>=H2){
            System.out.print("finalResult is maxCov = "+ H1);
            return resultHs1;
        }else {
            System.out.print("finalResult is maxRatio = "+ H2);
            return resultHs2;
        }
    }

    public static int simpleGreedySubFunction2(LinkedHashMap<String, Integer> resultHs, HashMap<Integer, HashSet<Integer>> subGraph, boolean maxRatio ) {
        // store the node and its coverage incremental
        PriorityQueue<relaxIndexNode> totalResult = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small

//        LinkedHashMap<String, Integer> resultHs = new LinkedHashMap<>();
        HashSet<String> resultlist = new HashSet<>();
        double totalPrice = 0;
        int maxCoverage = 0;

        for (int id: subGraph.keySet()){
            String idd = String.valueOf(id);
            if (datasetPrice.get(idd)<= budget)
                if (maxRatio){
                    totalResult.add(new relaxIndexNode(idd, datasetIdMapSignature.get(idd).size()/datasetPrice.get(idd), datasetIdMapSignature.get(idd).size()));
                }else{
                    totalResult.add(new relaxIndexNode(idd, datasetIdMapSignature.get(idd).size()));
                }
        }
        if (totalResult.size()>0){
            relaxIndexNode re = totalResult.poll();
            if (maxRatio){
                resultHs.put(re.resultId, (int)re.getUb());
                maxCoverage = (int)re.getUb();  // coverage
            }else {
                resultHs.put(re.resultId, (int)re.getLb());
                maxCoverage = (int) re.getLb();
            }

            totalPrice = datasetPrice.get(re.resultId);
            resultlist.add(re.resultId);
            int iteration = 1;
//            System.out.println("iteration "+iteration +" = "+re.resultId+","+ integerMaptoDatasetId.get(Integer.parseInt(re.resultId)) +", coverage = "+ maxCoverage +", getCurrentCost(resultPath) = "+totalPrice); //33.
            HashSet<Long> C = new HashSet<>();
            ArrayList<Long> arrayList = datasetIdMapSignature.get(re.resultId);
            for (int i=0; i< arrayList.size(); i++){
                C.add(arrayList.get(i));
            }

            boolean flag = true;
            while (flag && totalPrice <= budget) {
                iteration ++;
                PriorityQueue<relaxIndexNode> hs = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
                    @Override
                    public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                        return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
                    }
                }); // large -> small,   store datasetId and coverage incremental;
                for (String s : resultlist) {
                    int id = Integer.parseInt(s);
//                    System.out.println("s = "+s+", "+subGraph.containsKey(id));
                    if (subGraph.containsKey(id)) {
                        HashSet<Integer> edgeInfo = subGraph.get(id);
//                        System.out.println("edgeInfo = "+ edgeInfo);
                        for (int in : edgeInfo) {
                            String i = String.valueOf(in);
                            int setIntersection = datasetIdMapSignature.get(i).size();
                            if (!resultlist.contains(i) && (totalPrice+datasetPrice.get(i)) <= budget) {
                                for (long l: datasetIdMapSignature.get(i)){
                                    if (C.contains(l)){
                                        setIntersection--;
                                    }
                                }
                                if (setIntersection >0){
                                    if (maxRatio){
                                        hs.add(new relaxIndexNode(i, setIntersection/datasetPrice.get(i), setIntersection));
//                                        System.out.println(maxRatio+", "+i+", setIntersection = "+ setIntersection+", datasetPrice = "+datasetPrice.get(i));
                                    }else {
                                        hs.add(new relaxIndexNode(i, setIntersection));  //store id and coverage incremental
//                                        System.out.println(maxRatio+", "+ i+", setIntersection = "+ setIntersection+", datasetPrice = "+datasetPrice.get(i));
                                    }
                                }
                            }
                        }
                    }
                }
//                System.out.println("hs.size() = "+hs.size());

                if (hs.size()>0){
                    re = hs.peek();
                    if (maxRatio){
                        maxCoverage += re.ub;
                    }else {
                        maxCoverage += re.lb;
                    }
                    resultlist.add(re.resultId);
                    double price = datasetPrice.get(re.resultId) ;
                    totalPrice += price;
//                    if (maxRatio)
//                        System.out.println("iteration "+iteration +" = "+re.resultId+ "," + integerMaptoDatasetId.get(Integer.parseInt(re.resultId)) +", coverage = "+ re.ub +", getCurrentCost(resultPath) = "+totalPrice); //33.
//                    else
//                        System.out.println("iteration "+iteration +" = "+re.resultId+"," + integerMaptoDatasetId.get(Integer.parseInt(re.resultId)) +", coverage = "+ re.lb +", getCurrentCost(resultPath) = "+totalPrice); //33.
//                Double.parseDouble(String.format("%.2f", datasetPrice.get(re.resultId)));
                    for (int i = 0; i<datasetIdMapSignature.get(re.resultId).size(); i++){
                        C.add(datasetIdMapSignature.get(re.resultId).get(i));
                    }

                }else {
                    flag = false;
                }
            }
            if (maxRatio){
                System.out.println( "coverage = "+maxCoverage+ ",   the total price = "+ totalPrice);
            }else {
                System.out.println( "coverage = "+maxCoverage+ ",   the total price = "+ totalPrice);
            }
        }else {
            maxCoverage = 0;
        }


        return maxCoverage;

    }

    public static void getPminPmax(HashSet<String> subGraph){
        p_min = 1000000000000000000.0;
        p_max = 0;
        for (String s : subGraph){
            if (p_min > datasetPrice.get(s)){
                p_min = datasetPrice.get(s);
            }
            if (p_max < datasetPrice.get(s)){
                p_max = datasetPrice.get(s);
            }
        }
//        System.out.println("p_max = "+p_max+", "+"p_min = "+p_min);
    }

    public static relaxIndexNode findCenterRadius(HashMap<String, HashSet<String>> subGraph, double budget){
        PriorityQueue<relaxIndexNode> PQ = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                if (o1.getLb()> o2.getLb()) // small - large
                    return 1;
                else if (o1.getLb()< o2.getLb())
                    return -1;
                else {
                    if (o2.getUb() > o1.getUb()) // large -> small
                        return 1;
                    else
                        return -1;
                }
            }
        });


        for (String s: subGraph.keySet()){
            Queue<relaxIndexNode> queue = new LinkedList<relaxIndexNode>();
            HashSet<String> visitedNodeList = new HashSet<>();
            int depth = 0;
//            System.out.println("S = "+s);
            queue.add(new relaxIndexNode(s, depth));//
            visitedNodeList.add(s);
            while (!queue.isEmpty()){
                relaxIndexNode node = queue.poll();
                depth = (int) node.lb;
                HashSet<String> edgeInfo = subGraph.get(node.resultId);
                depth = depth+1;
                if (edgeInfo.size()>0){
                    for (String child: edgeInfo){
                        if (!visitedNodeList.contains(String.valueOf(child))){
                            queue.add(new relaxIndexNode(child, depth));
                            visitedNodeList.add(child);
                        }
                    }
                }

            }
            int nodeSize = visitedNodeList.size();
            PQ.add(new relaxIndexNode(s, depth-1, datasetIdMapSignature.get(s).size())); // radius from small->large, coverage from large->small
//            System.out.print("s= "+s + ", de= "+(depth-1)+", ");
        }
//        System.out.println();


        relaxIndexNode firstCenter = PQ.poll();
        relaxIndexNode re = new relaxIndexNode(firstCenter.resultId, firstCenter.getLb());
        String center = re.resultId;
        int radius = (int)firstCenter.getLb();

//        System.out.print("center*= "+center+",  ");
//        System.out.print("radius*= "+radius+",  ");
//        System.out.print("price= "+ outputPrice(datasetPrice.get(center))+",  ");
//        System.out.print("B need>= $"+ Math.ceil(radius*p_max+datasetPrice.get(center))+", ");

        while (!PQ.isEmpty()){
            if (datasetPrice.get(firstCenter.resultId)<= budget){
                break;
            }
            firstCenter = PQ.poll();
        }


        re = new relaxIndexNode(firstCenter.resultId, firstCenter.getLb());

//        System.out.println("center = "+ center+", "+ re.resultId);

//        if (re.resultId != center){
////            System.out.println();
//            System.out.print("selected center= "+re.resultId+",  ");
//            System.out.print("radius= "+(int)re.getLb()+",  ");
//            System.out.print("price= "+ String.format("%.1f",datasetPrice.get(re.resultId))+",  ");
//            System.out.println("B need>= $"+ Math.ceil(radius*p_max+datasetPrice.get(re.resultId))+", ");
//        }
//        System.out.println();

        return re;
    }
    public static relaxIndexNode fastFindCenterRadius(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph, double budget){
        PriorityQueue<relaxIndexNode> PQ = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small

        for (String s: subGraph.keySet()){
            Queue<relaxIndexNode> queue = new LinkedList<relaxIndexNode>();
            HashSet<String> visitedNodeList = new HashSet<>();
            int depth = 0;
            queue.add(new relaxIndexNode(s, depth));// 根节点，深度为1
            relaxIndexNode deepestNode = new relaxIndexNode(s, depth);
            visitedNodeList.add(s);
            while (!queue.isEmpty()){
                relaxIndexNode node = queue.poll();
                deepestNode = node;
                depth = (int) node.lb;
//                System.out.println(node.resultId+", "+depth);
                HashSet<String> edgeInfo = subGraph.get(node.resultId);
                depth = depth+1;
                if (edgeInfo.size()>0){
                    for (String child: edgeInfo){
                        if (!visitedNodeList.contains(child)){
                            queue.add(new relaxIndexNode(child, depth));
                            visitedNodeList.add(child);
                        }
                    }
                }
            }
            int nodeSize = visitedNodeList.size();
            PQ.add(new relaxIndexNode(deepestNode.resultId, depth-1, nodeSize));
            break;
        }
        relaxIndexNode deepestNode = PQ.poll();
        String deepthNodeId = deepestNode.resultId;
//        System.out.println("deepestNode = "+deepthNodeId);
        HashMap<String, indexNode> leafNodeList = constructBFS(subGraphIndexNodes, subGraph, deepthNodeId);
        int diameter = 0;
        for (String id: leafNodeList.keySet()){
            indexNode leafNode = leafNodeList.get(id);
            if (diameter < leafNode.getRadius() ){
                diameter = (int)leafNode.getRadius();
            }
        }
        int radius = (int)Math.floor((diameter+1)/2);
//        System.out.println("diameter = "+ diameter +", rrradius = "+radius);
        Queue<indexNode> queue = new LinkedList<indexNode>();
        HashMap<Integer, LinkedList<relaxIndexNode>> allLevelNodeList = new HashMap<>();
        queue.add(subGraphIndexNodes.get(deepthNodeId));
        while (!queue.isEmpty()){
            indexNode in = queue.poll();
            Set<indexNode> centerList = in.getNodelist();
            LinkedList<relaxIndexNode> levelNodeList = allLevelNodeList.getOrDefault((int)in.getRadius()+1, new LinkedList<>());
//            System.out.println("in.radius = "+in.getRadius()+",  in.deviceId ="+ in.getDeviceID()+", "+levelNodeList.size());
            for (indexNode i: centerList){
                levelNodeList.add(new relaxIndexNode(i.getDeviceID(), i.getRadius()));
                queue.add(i);
//                System.out.println(", add ="+i.getDeviceID()+", "+i.getRadius() +"; ");
            }
//            System.out.println("vvv = "+ in.getDeviceID()+", levelNodeList.size = "+levelNodeList.size() );
            allLevelNodeList.put((int)in.getRadius()+1, levelNodeList);
        }


        PriorityQueue<relaxIndexNode> centerPQ = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                if (o1.getLb()> o2.getLb()) // small - large
                    return 1;
                else if (o1.getLb()< o2.getLb())
                    return -1;
                else {
                    if (o2.getUb() > o1.getUb()) // large -> small
                        return 1;
                    else
                        return -1;
                }
            }
        });

        if (diameter % 2 == 0) {
            for (int r: allLevelNodeList.keySet()){
                LinkedList<relaxIndexNode> centerNodeList = allLevelNodeList.get(r);
                for (relaxIndexNode in: centerNodeList){
//                System.out.println("www = "+in.getResultId()+", "+ Math.floor(Math.abs(r-(diameter+1)/2))+", "+ datasetIdMapSignature.get(in.getResultId()).size());
                    centerPQ.add(new relaxIndexNode(in.getResultId(), Math.floor(Math.abs(r-(diameter+1)/2)), datasetIdMapSignature.get(in.getResultId()).size()));
                }
            }
        }else {
            for (int r: allLevelNodeList.keySet()){
                LinkedList<relaxIndexNode> centerNodeList = allLevelNodeList.get(r);
                for (relaxIndexNode in: centerNodeList){
//                System.out.println("www = "+in.getDeviceID()+", "+ Math.floor(Math.abs(r-(diameter+1)/2))+", "+ datasetIdMapSignature.get(in.getDeviceID()).size());
                    centerPQ.add(new relaxIndexNode(in.getResultId(), Math.floor(Math.abs(r-radius)), datasetIdMapSignature.get(in.getResultId()).size()));
                }
            }
        }

//        System.out.println("centerPQ.size() = "+centerPQ.size());
        relaxIndexNode centerIndexNode = centerPQ.poll();
        relaxIndexNode re = new relaxIndexNode(centerIndexNode.resultId, centerIndexNode.getLb()+radius);

        String center = re.resultId;
        radius = (int) re.getLb();
//        System.out.print("center*= "+center+",  ");
//        System.out.print("radius*= "+radius+",  ");
//        System.out.print("price= "+ outputPrice(datasetPrice.get(center))+",  ");
//        System.out.print("B need>= $"+ Math.ceil(radius*p_max+datasetPrice.get(center))+", ");

        while (!centerPQ.isEmpty()){
            if (datasetPrice.get(centerIndexNode.resultId)<= budget){
                break;
            }
            centerIndexNode = centerPQ.poll();
        }
        re = new relaxIndexNode(centerIndexNode.resultId, centerIndexNode.getLb()+radius);
//        if (re.resultId != center){
//            System.out.println();
//            System.out.print("selected center= "+re.resultId+",  ");
//            System.out.print("radius= "+(int)re.getLb()+",  ");
//            System.out.print("price= "+ String.format("%.1f",datasetPrice.get(re.resultId))+",  ");
//            System.out.println("B need>= $"+ Math.ceil(radius*p_max+datasetPrice.get(re.resultId))+", ");
//        }

//        System.out.println();

        return re;
    }


    public static HashMap<String, indexNode> constructBFS(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph, String id) {
//        1,construct the BFS tree
        for (String s:subGraphIndexNodes.keySet()){
            subGraphIndexNodes.get(s).setDeviceID(s);
            subGraphIndexNodes.get(s).setroot(-2); //lable it as internal node
            subGraphIndexNodes.get(s).setRadius(0);
            subGraphIndexNodes.get(s).getNodelist().clear();
        }

        double radius = 0.0;
        BFS_root = subGraphIndexNodes.get(id);
        BFS_root.setRadius(radius);
        Queue<indexNode> queue = new LinkedList<indexNode>();
        HashSet<String> visitedNodeList = new HashSet<>();
        HashMap<String, indexNode> leafNodeList = new HashMap<>();
        queue.add(BFS_root);//
        visitedNodeList.add(id);
        while (!queue.isEmpty()){
            indexNode parant = queue.poll();
            HashSet<String> edgeInfo = subGraph.get(parant.getDeviceID());
//            System.out.println("parant = "+ parant.deviceID+ ", depth = "+parant.getRadius()+", edgeInfo.size() = "+edgeInfo.size());
            if (edgeInfo.size()>0){
                boolean isLeaf = true;
                for (String child: edgeInfo){
                    if (!visitedNodeList.contains(child)){
                        subGraphIndexNodes.get(child).setParentNode(parant);
                        subGraphIndexNodes.get(child).setRadius(parant.getRadius()+1);
//                        System.out.println("child = "+child+", depth = "+datalakeRootNodes.get(child).getRadius());
                        parant.addNodes(subGraphIndexNodes.get(child));
                        queue.add(subGraphIndexNodes.get(child));
                        visitedNodeList.add(child);
                        isLeaf = false;
//                        parant.setroot(-2);
                    }
                }
                if (isLeaf){
//                        System.out.println("1. "+ parant.getDeviceID()+" is set -1");
                    parant.setroot(-1); //lable it as leaf node
                }
            }else {
//                System.out.println("2. "+parant.getDeviceID()+" is set -1");
                parant.setroot(-1);
            }
        }

        for (String s:subGraphIndexNodes.keySet()){
            if (subGraphIndexNodes.get(s).getroot() == -1){
                leafNodeList.put(s, subGraphIndexNodes.get(s));
//                System.out.println("leafNodeList.add "+ s);
            }
        }
//        System.out.println("leafNodeList.size() = "+leafNodeList.size());
        return leafNodeList;
    }

    public static double afterAddNewNodesPrice(HashSet<String> resultPath, HashMap<String, ArrayList<String>> recordPaths, String id){
        ArrayList<String>  path = recordPaths.get(id);
        double currentCost = 0;
        for (String s: path){
            if (!resultPath.contains(s)){
                currentCost += datasetPrice.get(s);
            }
        }
        return currentCost;
    }

    public static double afterAddNewNodesPrice(HashSet<String> resultPath,  ArrayList<String> path){
        double currentCost = 0;
        for (String s: path){
            if (!resultPath.contains(s)){
                currentCost += datasetPrice.get(s);
            }
        }
        return currentCost;
    }

    public static double afterAddNewNodesPrice(HashSet<String> resultPath, String node){
        double currentCost = 0;

        if (!resultPath.contains(node)){
            currentCost += datasetPrice.get(node);
        }

        return currentCost;
    }

    public static double getCurrentCost(HashSet<String> resultPath){
        double currentCost = 0;
        if (resultPath.size()>0){
            for (String s : resultPath){
                currentCost += datasetPrice.get(s);
            }
        }
        return currentCost;
    }


    public static int findOneConnectedGraph(HashMap<Integer, HashSet<Integer>> graphInfo){

        HashSet<String> resultPath = new HashSet<>();//记录结果集合的ID
        HashSet<Long> C = new HashSet<>(); //记录结果集合覆盖的元素ID
        PriorityQueue<relaxIndexNode> totalResult = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small
        for (int id: graphInfo.keySet()){
            String idd = String.valueOf(id);
            if (datasetPrice.get(idd)<= budget){
                totalResult.add(new relaxIndexNode(idd, datasetIdMapSignature.get(idd).size()));
            }
        }

        if (totalResult.size()>0) {
            relaxIndexNode centerIndexNode = totalResult.poll();
            while (!totalResult.isEmpty()){
                if (datasetPrice.get(centerIndexNode.resultId)<= budget){
                    break;
                }
                centerIndexNode = totalResult.poll();
            }
            String selctedStringId = centerIndexNode.resultId;
            if (datasetPrice.get(selctedStringId)<=budget){
                resultPath.add(selctedStringId); // add root node;
                C.addAll(datasetIdMapSignature.get(selctedStringId));
                boolean flag = true;
                while (getCurrentCost(resultPath)<budget && flag) {
                    // random choose a connected node
                    HashSet<Integer> allConnectedNodes = new HashSet<>();
                    int i = Integer.parseInt(selctedStringId);
                    HashSet<Integer> connectedNodes = graphInfo.get(i);
                    allConnectedNodes.addAll(connectedNodes);

                    String tempStringId = selctedStringId;
                    for (int selectedIntegerId: allConnectedNodes){
                        String stringId = String.valueOf(selectedIntegerId);
                        if (getCurrentCost(resultPath)+afterAddNewNodesPrice(resultPath,stringId)<budget && !resultPath.contains(stringId)){
                            selctedStringId = stringId;
                            tempStringId = selctedStringId;
                            break;
                        }
                    }
                    if (selctedStringId.equals(tempStringId)){
                        flag = false;
                    }else {
                        resultPath.add(selctedStringId);
                        allConnectedNodes.remove(Integer.parseInt(selctedStringId));
                    }
                }
            }
        }
        return C.size();
    }


    public static int findOneConnectedGraphForallNodes(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> graphInfo){

        HashSet<String> resultPath = new HashSet<>();//记录结果集合的ID
        HashSet<Long> C = new HashSet<>(); //记录结果集合覆盖的元素ID

        PriorityQueue<relaxIndexNode> totalResult = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small

        ArrayList<String> nodeList = new ArrayList<>();
        for (String idd: graphInfo.keySet()){
            if (datasetPrice.get(idd)<= budget){
                nodeList.add(idd);
            }
        }
        for (String selctedStringId: nodeList) {
//            int number = (int) (Math.random() * (nodeList.size()));
//            String selctedStringId = nodeList.get(number);
            C.clear();
            resultPath.clear();
            resultPath.add(selctedStringId); // add root node;
            System.out.print("first = "+selctedStringId+","+integerMaptoDatasetId.get(Integer.parseInt(selctedStringId))+"; ");
            C.addAll(datasetIdMapSignature.get(selctedStringId));
            HashSet<String> allConnectedNodes = new HashSet<>();
            boolean flag = true;
            while (getCurrentCost(resultPath) < budget && flag) {
                HashSet<String> connectedNodes = graphInfo.get(selctedStringId);
                allConnectedNodes.addAll(connectedNodes);
                String tempStringId = selctedStringId;
                for (String stringId: allConnectedNodes){
                    if (getCurrentCost(resultPath)+afterAddNewNodesPrice(resultPath,stringId)<budget && !resultPath.contains(stringId)){
                        selctedStringId = stringId;
//                        tempStringId = selctedStringId;
                        break;
                    }
                }
                if (selctedStringId.equals(tempStringId)){
                    flag = false;
                }else {
                    resultPath.add(selctedStringId);
//                    System.out.println("resultPath = "+ resultPath);
                    allConnectedNodes.remove(selctedStringId);
                    System.out.print(selctedStringId+","+integerMaptoDatasetId.get(Integer.parseInt(selctedStringId))+"; ");
                }
            }
//            System.out.println("C.size = "+C.size());
            totalResult.add(new relaxIndexNode(selctedStringId, C.size()));
            break;
        }

        int coverage = 0;
        if (totalResult.size()>0){
            coverage = (int)totalResult.peek().lb;
        }
        return coverage;
    }


    public static int findOneConnectedGraphRandomNode(HashMap<String, indexNode> subGraphIndexNodes, HashMap<Integer, HashSet<Integer>> graphInfo){

        HashSet<String> resultPath = new HashSet<>();//记录结果集合的ID
        HashSet<Long> C = new HashSet<>(); //记录结果集合覆盖的元素ID
        int coverage = 0;
        PriorityQueue<relaxIndexNode> totalResult = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small

        ArrayList<String> nodeList = new ArrayList<>();
        for (int id: graphInfo.keySet()){
            String idd = String.valueOf(id);
            if (datasetPrice.get(idd)<= budget){
                nodeList.add(idd);
            }
        }

        int number = (int) (Math.random() * (nodeList.size()));
        String selctedStringId = nodeList.get(number);
        C.clear();
        resultPath.clear();
        resultPath.add(selctedStringId); // add root node;
        C.addAll(datasetIdMapSignature.get(selctedStringId));
        boolean flag = true;
        while (getCurrentCost(resultPath) < budget && flag) {
            // random choose a connected node
            HashSet<Integer> allConnectedNodes = new HashSet<>();
            int i = Integer.parseInt(selctedStringId);
            HashSet<Integer> connectedNodes = graphInfo.get(i);
            allConnectedNodes.addAll(connectedNodes);
//                System.out.println("allConnectedNodes = "+allConnectedNodes);
            String tempStringId = selctedStringId;
            for (int selectedIntegerId: allConnectedNodes){
                String stringId = String.valueOf(selectedIntegerId);
                if (getCurrentCost(resultPath)+afterAddNewNodesPrice(resultPath,stringId)<budget && !resultPath.contains(stringId)){
                    selctedStringId = stringId;
                    tempStringId = selctedStringId;
                    break;
                }
            }
            if (selctedStringId.equals(tempStringId)){
                flag = false;
            }else {
                resultPath.add(selctedStringId);
//                    System.out.println("resultPath = "+ resultPath);
                allConnectedNodes.remove(Integer.parseInt(selctedStringId));
            }
//                System.out.println("selctedStringId = "+selctedStringId+", lastSelectID = "+tempStringId);
        }
//            System.out.println("C.size = "+C.size());


        return C.size();
    }


    public static int CombinatorialAlgorithm(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph) {
        HashSet<String> resultPath = new HashSet<>();//记录结果集合的ID
        HashSet<Long> C = new HashSet<>(); //记录结果集合覆盖的元素ID
        HashSet<String> tempResultPath = new HashSet<>(); //记录每次迭代的结果，所以每次迭代之前需要clear
        HashSet<Long> tempC = new HashSet<>();

        for (String center: subGraphIndexNodes.keySet()) {
//            System.out.println("v = "+center);
            HashMap<String, indexNode> leafNodeList = constructBFS(subGraphIndexNodes, subGraph, center);
            tempResultPath.clear(); //                remember clear
            tempC.clear();
            if (datasetPrice.get(center)<budget) {
                //            2,记录S_L, Path,和 Price
                HashMap<String, HashSet<Long>> S_L = new HashMap<>();
                HashMap<String, ArrayList<String>> recordPaths = new HashMap<>(); // leafnodeId, <leaf, leaf.parent, leaf.parent.parent, ...>
                HashMap<String, Double> recordPrice = new HashMap<>();
                for (String leafId : leafNodeList.keySet()) {
                    indexNode leaf = leafNodeList.get(leafId);
                    HashSet<Long> hs = new HashSet<>();
                    ArrayList<String> path = new ArrayList<>();
                    double S_L_price = datasetPrice.get(leaf.getDeviceID());
                    ArrayList<Long> arrayList = datasetIdMapSignature.get(leaf.getDeviceID());
                    path.add(leaf.getDeviceID());
                    for (int i = 0; i < arrayList.size(); i++)
                        hs.add(arrayList.get(i));

                    indexNode parant = leaf.getParentNode();
                    while (!parant.deviceID.equals(center)) {
                        S_L_price += datasetPrice.get(parant.getDeviceID());
                        arrayList = datasetIdMapSignature.get(parant.getDeviceID());
                        path.add(parant.getDeviceID());
                        for (int i = 0; i < arrayList.size(); i++) {
                            hs.add(arrayList.get(i));
                        }
                        parant = parant.getParentNode();
                    }
                    S_L.put(leaf.getDeviceID(), hs);
                    recordPaths.put(leaf.getDeviceID(), path);
                    recordPrice.put(leaf.getDeviceID(), S_L_price);
                }

//                System.out.println("* = "+center);
//                for (String leafId : leafNodeList.keySet()){
//                    System.out.println("leafId = "+leafId);
//                }

                tempResultPath.add(center); // add root node;
                for (long l: datasetIdMapSignature.get(center))
                    tempC.add(l);
                boolean flag = true;
                while (getCurrentCost(tempResultPath)<budget && flag) {
                    double currentBestratio = 0.0;
                    ArrayList<String> currentBestPath = new ArrayList<>();

                    for (String leafnodeId : recordPaths.keySet()) {
                        ArrayList<String> leaftorootPath = recordPaths.get(leafnodeId);
                        // for each u, find the path with maximum ratio;
                        for (int i = leaftorootPath.size() - 1; i >= 0; i--) {
                            String u = leaftorootPath.get(i);
                            if (!tempResultPath.contains(u)) {
                                ArrayList<String> tempPath = new ArrayList<>();
                                for (int j = leaftorootPath.size() - 1; j >= i; j--)
                                    tempPath.add(leaftorootPath.get(j));
                                if (getCurrentCost(tempResultPath) + afterAddNewNodesPrice(tempResultPath, tempPath) <= budget) {
                                    HashSet<Long> P_vu = new HashSet<>();
                                    for (String s : tempPath) {
                                        for (long l : datasetIdMapSignature.get(s)) {
                                            if (!tempC.contains(l))
                                                P_vu.add(l);
                                        }
                                    }
                                    HashSet<String> l_vu = new HashSet<>();
                                    for (String s : tempPath) {
                                        if (!tempResultPath.contains(s))
                                            l_vu.add(s);
                                    }

                                    double tempThresh = 0.0;
                                    if (l_vu.size() > 0)
                                        tempThresh = P_vu.size() / l_vu.size();  //计算

                                    if (currentBestratio < tempThresh) {
                                        currentBestPath.clear();
                                        for (String s : tempPath) {
                                            currentBestPath.add(s);
                                        }
                                        currentBestratio = tempThresh;
                                    }

//                                    System.out.println("u = " + u + ", "+ "currentBestratio = " + currentBestratio + ", " + " tempThresh = " + tempThresh + ", "+ "currentBestPath = " + currentBestPath);

                                }
                            }
                        }
                    }
                    if (currentBestratio == 0) {
                        flag = false;
//                            System.out.println("currentBestratio==0");
                    } else {
                        for (String node : currentBestPath) {
                            tempResultPath.add(node);
//                                System.out.println("nodeeee = "+node);
                            for (long l : datasetIdMapSignature.get(node)) {
                                if (!tempC.contains(l))
                                    tempC.add(l);
                            }
                        }
                    }
                }
            }

        }
        if (tempC.size()> C.size()){
            for (String s: tempResultPath){
                resultPath.add(s);
            }
            // note address reference;
//                C.clear();
            for (long l: tempC)
                C.add(l);
        }



//        System.out.println("resultPath = ");
//        for (String s: resultPath){
//            System.out.print(s+", ");
//        }
//        System.out.print("budget = "+getCurrentCost(resultPath)+", ");

        return C.size();

    }

    public static int CombinatorialAlgorithmSpeedMaxGain(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph) {
        HashSet<String> resultPath = new HashSet<>();//记录结果集合的ID
        HashSet<Long> C = new HashSet<>(); //记录结果集合覆盖的元素ID
        HashSet<String> tempResultPath = new HashSet<>(); //记录每次迭代的结果，所以每次迭代之前需要clear
        HashSet<Long> tempC = new HashSet<>();

        PriorityQueue<relaxIndexNode> totalResult = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
            @Override
            public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
            }
        }); // large -> small
        for (String idd: subGraphIndexNodes.keySet()){
            if (datasetPrice.get(idd)<= budget){
                totalResult.add(new relaxIndexNode(idd, datasetIdMapSignature.get(idd).size()));
            }
        }

        if (totalResult.size()>0) {
            relaxIndexNode centerIndexNode = totalResult.poll();
            while (!totalResult.isEmpty()){
                if (datasetPrice.get(centerIndexNode.resultId)<= budget){
                    break;
                }
                centerIndexNode = totalResult.poll();
            }
            String center = centerIndexNode.resultId;
//            System.out.println("v = "+center);
            HashMap<String, indexNode> leafNodeList = constructBFS(subGraphIndexNodes, subGraph, center);
            tempResultPath.clear(); //                remember clear
            tempC.clear();
            if (datasetPrice.get(center)<budget) {
                //            2,记录S_L, Path,和 Price
                HashMap<String, HashSet<Long>> S_L = new HashMap<>();
                HashMap<String, ArrayList<String>> recordPaths = new HashMap<>(); // leafnodeId, <leaf, leaf.parent, leaf.parent.parent, ...>
                HashMap<String, Double> recordPrice = new HashMap<>();
                for (String leafId : leafNodeList.keySet()) {
                    indexNode leaf = leafNodeList.get(leafId);
                    HashSet<Long> hs = new HashSet<>();
                    ArrayList<String> path = new ArrayList<>();
                    double S_L_price = datasetPrice.get(leaf.getDeviceID());
                    ArrayList<Long> arrayList = datasetIdMapSignature.get(leaf.getDeviceID());
                    path.add(leaf.getDeviceID());
                    for (int i = 0; i < arrayList.size(); i++)
                        hs.add(arrayList.get(i));

                    indexNode parant = leaf.getParentNode();
                    while (!parant.deviceID.equals(center)) {
                        S_L_price += datasetPrice.get(parant.getDeviceID());
                        arrayList = datasetIdMapSignature.get(parant.getDeviceID());
                        path.add(parant.getDeviceID());
                        for (int i = 0; i < arrayList.size(); i++) {
                            hs.add(arrayList.get(i));
                        }
                        parant = parant.getParentNode();
                    }
                    S_L.put(leaf.getDeviceID(), hs);
                    recordPaths.put(leaf.getDeviceID(), path);
                    recordPrice.put(leaf.getDeviceID(), S_L_price);
                }

//                System.out.println("* = "+center);
//                for (String leafId : leafNodeList.keySet()){
//                    System.out.println("leafId = "+leafId);
//                }

                tempResultPath.add(center); // add root node;
                for (long l: datasetIdMapSignature.get(center))
                    tempC.add(l);
                boolean flag = true;
                while (getCurrentCost(tempResultPath)<budget && flag) {
                    double currentBestratio = 0.0;
                    ArrayList<String> currentBestPath = new ArrayList<>();

                    for (String leafnodeId : recordPaths.keySet()) {
                        ArrayList<String> leaftorootPath = recordPaths.get(leafnodeId);
                        // for each u, find the path with maximum ratio;
                        for (int i = leaftorootPath.size() - 1; i >= 0; i--) {
                            String u = leaftorootPath.get(i);
                            if (!tempResultPath.contains(u)) {
                                ArrayList<String> tempPath = new ArrayList<>();
                                for (int j = leaftorootPath.size() - 1; j >= i; j--)
                                    tempPath.add(leaftorootPath.get(j));
                                if (getCurrentCost(tempResultPath) + afterAddNewNodesPrice(tempResultPath, tempPath) <= budget) {
                                    HashSet<Long> P_vu = new HashSet<>();
                                    for (String s : tempPath) {
                                        for (long l : datasetIdMapSignature.get(s)) {
                                            if (!tempC.contains(l))
                                                P_vu.add(l);
                                        }
                                    }
                                    HashSet<String> l_vu = new HashSet<>();
                                    for (String s : tempPath) {
                                        if (!tempResultPath.contains(s))
                                            l_vu.add(s);
                                    }

                                    double tempThresh = 0.0;
                                    if (l_vu.size() > 0)
                                        tempThresh = P_vu.size() / l_vu.size();  //计算

                                    if (currentBestratio < tempThresh) {
                                        currentBestPath.clear();
                                        for (String s : tempPath) {
                                            currentBestPath.add(s);
                                        }
                                        currentBestratio = tempThresh;
                                    }

//                                    System.out.println("u = " + u + ", "+ "currentBestratio = " + currentBestratio + ", " + " tempThresh = " + tempThresh + ", "+ "currentBestPath = " + currentBestPath);

                                }
                            }
                        }
                    }
                    if (currentBestratio == 0) {
                        flag = false;
//                            System.out.println("currentBestratio==0");
                    } else {
                        for (String node : currentBestPath) {
                            tempResultPath.add(node);
//                                System.out.println("nodeeee = "+node);
                            for (long l : datasetIdMapSignature.get(node)) {
                                if (!tempC.contains(l))
                                    tempC.add(l);
                            }
                        }
                    }
                }
            }

        }
        if (tempC.size()> C.size()){
            for (String s: tempResultPath){
                resultPath.add(s);
            }
            // note address reference;
//                C.clear();
            for (long l: tempC)
                C.add(l);
        }


//        System.out.print("C.size = "+ C.size()+": ");
//        for (String s: resultPath){
//            System.out.print(s+", "+integerMaptoDatasetId.get(Integer.parseInt(s))+"; ");
//        }
//        System.out.println();
//        System.out.print("budget = "+getCurrentCost(resultPath)+", ");
        return C.size();

    }

    public static int CombinatorialAlgorithmSpeedMaxSize(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph) {
        HashSet<String> resultPath = new HashSet<>();
        HashSet<Long> C = new HashSet<>();
        HashSet<String> tempResultPath = new HashSet<>();
        HashSet<Long> tempC = new HashSet<>();

        HashSet<String> totalResult = new HashSet<>();
        for (String idd: subGraphIndexNodes.keySet()){
            if (datasetPrice.get(idd)<= budget){
                totalResult.add(idd);
            }
        }
//        System.out.println("totalResult = "+totalResult );

        if (totalResult.size()>0) {
            String center = "";
            for (String s: totalResult){
                center = s;
                break;
            }
//            System.out.println("v = "+center);
            HashMap<String, indexNode> leafNodeList = constructBFS(subGraphIndexNodes, subGraph, center);
            tempResultPath.clear(); //                remember clear
            tempC.clear();
            if (datasetPrice.get(center)<budget) {
                //            2,记录S_L, Path,和 Price
                HashMap<String, HashSet<Long>> S_L = new HashMap<>();
                HashMap<String, ArrayList<String>> recordPaths = new HashMap<>(); // leafnodeId, <leaf, leaf.parent, leaf.parent.parent, ...>
                HashMap<String, Double> recordPrice = new HashMap<>();
                for (String leafId : leafNodeList.keySet()) {
                    indexNode leaf = leafNodeList.get(leafId);
                    HashSet<Long> hs = new HashSet<>();
                    ArrayList<String> path = new ArrayList<>();
                    double S_L_price = datasetPrice.get(leaf.getDeviceID());
                    ArrayList<Long> arrayList = datasetIdMapSignature.get(leaf.getDeviceID());
                    path.add(leaf.getDeviceID());
                    for (int i = 0; i < arrayList.size(); i++)
                        hs.add(arrayList.get(i));

                    indexNode parant = leaf.getParentNode();
                    while (!parant.deviceID.equals(center)) {
                        S_L_price += datasetPrice.get(parant.getDeviceID());
                        arrayList = datasetIdMapSignature.get(parant.getDeviceID());
                        path.add(parant.getDeviceID());
                        for (int i = 0; i < arrayList.size(); i++) {
                            hs.add(arrayList.get(i));
                        }
                        parant = parant.getParentNode();
                    }
                    S_L.put(leaf.getDeviceID(), hs);
                    recordPaths.put(leaf.getDeviceID(), path);
                    recordPrice.put(leaf.getDeviceID(), S_L_price);
                }

                tempResultPath.add(center); // add root node;
                for (long l: datasetIdMapSignature.get(center))
                    tempC.add(l);
                boolean flag = true;
                int count = 0;
                while (getCurrentCost(tempResultPath)<budget && flag) {
                    count++;
//                    System.out.println("count = "+count);
                    double currentBestratio = 0.0;
                    ArrayList<String> currentBestPath = new ArrayList<>();
                    String maxCovPathId = "";
                    if (recordPaths.size()>0){
                        for (String leafnodeId : recordPaths.keySet()) {
                            ArrayList<String> leaftorootPath = recordPaths.get(leafnodeId);
                            if (getCurrentCost(tempResultPath) + afterAddNewNodesPrice(tempResultPath, leaftorootPath) <= budget) {
                                double tempThresh = 0;
                                for (String s : leaftorootPath) {
                                    tempThresh += datasetIdMapSignature.get(s).size();
                                }
                                if (currentBestratio < tempThresh) {
                                    maxCovPathId = leafnodeId;
                                    currentBestratio = tempThresh;
                                }
                            }
                        }
                    }

                    if (currentBestratio == 0 || maxCovPathId.equals("")) {
                        flag = false;
                    } else {
                        currentBestPath.clear();
//                        System.out.println("maxCovPathId = "+maxCovPathId);
//                        System.out.println("recordPaths = "+recordPaths);
                        for (String s : recordPaths.get(maxCovPathId)) {
                            currentBestPath.add(s);
                        }
                        for (String node : currentBestPath) {
                            tempResultPath.add(node);
                            for (long l : datasetIdMapSignature.get(node)) {
                                if (!tempC.contains(l))
                                    tempC.add(l);
                            }
                        }
                        recordPaths.remove(maxCovPathId);
                    }
                }
            }

        }
        if (tempC.size()> C.size()){
            for (String s: tempResultPath){
                resultPath.add(s);
            }
            // note address reference;
//                C.clear();
            for (long l: tempC)
                C.add(l);
        }


//        System.out.print("C.size = "+ C.size()+": ");
//        for (String s: resultPath){
//            System.out.print(s+", "+integerMaptoDatasetId.get(Integer.parseInt(s))+"; ");
//        }
//        System.out.println();

        return C.size();

    }


    public static LinkedHashMap<String, Integer> IBGPA(HashMap<String, indexNode> subGraphIndexNodes, HashMap<String, HashSet<String>> subGraph, boolean speedTrick, ArrayList<Long> findCenterTime) throws CloneNotSupportedException{
//        1,为每一个节点构建BFS树，然后找到center node和radius
//        注意的是如果图不是一个连通图，就要注意找到多个图
//         1,记录节点名字和层数，然后返回最小的层数，按照从小到大排序，
        HashMap<String, indexNode> leafNodeList = new HashMap<>();
        String center = "";
        int radius = 0;

        if (speedTrick == true){
            // 快速找center node
            long time03 = System.currentTimeMillis();
            relaxIndexNode centerRadius = fastFindCenterRadius(subGraphIndexNodes, subGraph, budget);
            long time04 = System.currentTimeMillis();
//            System.out.print("fast findCenter ="+(time04-time03) +"ms, ");
            findCenterTime.add((time04-time03));
            center = centerRadius.resultId;
            radius = (int)centerRadius.getLb();
            leafNodeList = constructBFS(subGraphIndexNodes, subGraph, center);
        }else {
            // 遍历找center node
            long time01 = System.currentTimeMillis();
            relaxIndexNode centerRadius = findCenterRadius(subGraph, budget);
            long time02 = System.currentTimeMillis();
//            System.out.print("find center ="+(time02-time01) +"ms, ");
            findCenterTime.add((time02-time01));
            center = centerRadius.resultId;
            radius = (int)centerRadius.getLb();
            leafNodeList = constructBFS(subGraphIndexNodes, subGraph, center); //            1,找到所有的叶子节点
        }


        LinkedHashMap<String, Integer> resultHs = new LinkedHashMap<>();
        HashSet<String> resultPath = new HashSet<>();

        if (datasetPrice.get(center)<budget){
            //            2,记录S_L, Path,和 Price
            HashMap<String, HashSet<Long>> S_L = new HashMap<>();
            HashMap<String, ArrayList<String>> recordPaths = new HashMap<>();
            HashMap<String, Double> recordPrice = new HashMap<>();
            for (String leafId: leafNodeList.keySet()){
                indexNode leaf = leafNodeList.get(leafId);
                HashSet<Long> hs = new HashSet<>();
                ArrayList<String> path = new ArrayList<>();
                double S_L_price = datasetPrice.get(leaf.getDeviceID());
                ArrayList<Long> arrayList = datasetIdMapSignature.get(leaf.getDeviceID());
                path.add(leaf.getDeviceID());
                for (int i = 0; i<arrayList.size(); i++){
                    hs.add(arrayList.get(i));
                }

                indexNode parant = leaf.getParentNode();
//                System.out.println("parant.deviceID = "+parant.deviceID);
                while (!parant.deviceID.equals(center)){
                    S_L_price += datasetPrice.get(parant.getDeviceID());
                    arrayList = datasetIdMapSignature.get(parant.getDeviceID());
                    path.add(parant.getDeviceID());
                    for (int i = 0; i<arrayList.size(); i++){
                        hs.add(arrayList.get(i));
                    }
                    parant = parant.getParentNode();
                }
                S_L.put(leaf.getDeviceID(), hs);
                recordPaths.put(leaf.getDeviceID(), path);
                recordPrice.put(leaf.getDeviceID(), S_L_price);
            }
//        System.out.println("recordPaths.size() = "+recordPaths.size());
//        for (String s: recordPaths.keySet()){
//            System.out.print("leaf node = "+ s+", price = " + recordPrice.get(s) +" : ");
//            ArrayList<String> ay = recordPaths.get(s);
//            for (String s2: ay){
//                System.out.print(s2+", "+ datasetPrice.get(s2)+"; ");
//            }
//            System.out.println();
//        }


//        System.out.println();
//        System.out.println("maximum coverage:");
            long ltime1 = System.currentTimeMillis();
            LinkedHashMap<String, Integer> resultHs1 = new LinkedHashMap<>();
            HashSet<Long> H_1 = UsedUpBudget(resultHs1, center, leafNodeList, S_L, recordPaths, recordPrice, false);
            long ltime2 = System.currentTimeMillis();
//            System.out.println("H_1.size() = "+H_1.size()+", time = "+ (ltime2-ltime1)+"ms ");
//        System.out.println("maximum coverage = "+H_1.size());

//        System.out.println("maximum ratio:");
            LinkedHashMap<String, Integer> resultHs2 = new LinkedHashMap<>();
            HashSet<Long> H_2 = UsedUpBudget(resultHs2, center, leafNodeList, S_L, recordPaths, recordPrice, true);
            long ltime3 = System.currentTimeMillis();
//            System.out.println("H_2.size() = "+H_2.size()+", time = "+ (ltime3-ltime2)+"ms ");
//        System.out.println("maximum ratio = "+H_2.size());
            if (H_1.size()> H_2.size()){
                for (String s: resultHs1.keySet())
                {
                    if (recordPaths.containsKey(s)){
                        ArrayList<String> ay = recordPaths.get(s);
                        for (String s2: ay)
                            resultPath.add(s2);
                    }else {
                        resultPath.add(s);
                    }
                }
//            System.out.print("remaingBudget = "+ outputPrice(budget - getCurrentCost(resultPath))+", ");
//            System.out.println("finalResult is maxCov = "+ H_1.size());
                resultHs = resultHs1;
//            for (String s: resultHs.keySet()){
//                System.out.println("s = "+s+",  coverage = "+resultHs.get(s));
//            }
            }else{
                for (String s: resultHs2.keySet())
                {
                    if (recordPaths.containsKey(s)){
                        ArrayList<String> ay = recordPaths.get(s);
                        for (String s2: ay)
                            resultPath.add(s2);
                    }else {
                        resultPath.add(s);
                    }
                }
//            System.out.print("remaingBudget = "+ outputPrice(budget - getCurrentCost(resultPath))+", ");
//            System.out.println("finalResult is maxRatio = "+ H_2.size());
                resultHs = resultHs2;
//            for (String s: resultHs.keySet()){
//                System.out.println("s = "+s+",  coverage = "+resultHs.get(s));
//            }
            }

        }else {
            resultHs.put(center, -1);
            resultPath.add(center);
        }

//        System.out.println("remaing budget = "+ outputPrice(budget - getCurrentCost(resultPath)));
//        if((budget - getCurrentCost(resultPath))<0){
//            System.out.print("budget overflow: ");
//            for (String s: resultPath){
//                System.out.print(s+",  ");
//            }
//            System.out.println();
//        }

        return resultHs;
    }

    public static HashSet<Long> UsedUpBudget(LinkedHashMap<String, Integer> rootleafNodeIncreHm, String v, HashMap<String, indexNode> leafnodelist, HashMap<String, HashSet<Long>> S_L, HashMap<String, ArrayList<String>> recordPaths, HashMap<String, Double> recordPrice, boolean maxRatio){

        String center = v;
        HashMap<String, indexNode> leafNodeList = new HashMap<>();
        for (String id: leafnodelist.keySet()){
            leafNodeList.put(id, leafnodelist.get(id));
        }
//            3，写while循环，找贪心解
        HashSet<String> resultPath = new HashSet<>();
        resultPath.add(center);
        rootleafNodeIncreHm.put(center, datasetIdMapSignature.get(center).size());
        ArrayList<Long> arrayList = datasetIdMapSignature.get(center);

        HashSet<Long> C = new HashSet<>();
        for (int i =0; i<arrayList.size(); i++){
            C.add(arrayList.get(i));
        }

        int iteration = 1;

//        System.out.println("iteration# "+iteration +" = "+center +", root node coverage = "+ datasetIdMapSignature.get(center).size() +", cost = "+datasetPrice.get(center));

        while (getCurrentCost(resultPath) < budget && leafNodeList.size()>0){  // &&  !leafNodeList.isEmpty()
//            选择最大max
            iteration++;
            PriorityQueue<relaxIndexNode> coverageIncremental = new PriorityQueue<>(new Comparator<relaxIndexNode>() {
                @Override
                public int compare(relaxIndexNode o1, relaxIndexNode o2) {
                    return (o2.getLb() - o1.getLb() > 0) ? 1 : -1;
                }
            }); // large -> small

            HashSet<String> leafDatasetId = new HashSet<>();
            for (String leafId: leafNodeList.keySet())
                leafDatasetId.add(leafId);
            for (String leafId: leafDatasetId){
                //                System.out.println("leafId = "+leafId);
                indexNode leaf = leafNodeList.get(leafId);
                HashSet<Long> S_L_i = S_L.get(leaf.getDeviceID());
                int setIntersection = S_L_i.size();
                for (Long l: S_L_i){
                    if (C.contains(l)){
                        setIntersection --;
                    }
                }


                if (setIntersection>0){
                    if (maxRatio == false){
                        coverageIncremental.add(new relaxIndexNode(leaf.getDeviceID(), setIntersection));
                    }else {
//                        1, 计算路径的价格，
                        ArrayList<String> path = recordPaths.get(leaf.getDeviceID());
                        HashSet<String> newNode = new HashSet<>();
                        for (int i=0; i<path.size(); i++){
                            String id = path.get(i);
                            if (!resultPath.contains(id)){
                                newNode.add(id);
                            }
                        }
                        coverageIncremental.add(new relaxIndexNode(leaf.getDeviceID(), (double) setIntersection /getCurrentCost(newNode) , setIntersection));
                    }
                }
            }// 按照leaf node的增量从大到小排序


//            System.out.println("maxRatio!!! = "+maxRatio);
//            for (relaxIndexNode re: coverageIncremental){
//                System.out.print(re.getResultId()+", "+ re.getLb()+" , "+ re.getUb()+";  ");
//            }
//            System.out.println();


//                控制选择的不能大于budget
//            System.out.println("coverageIncremental.size() = "+ coverageIncremental.size());
            if (coverageIncremental.size()>0){
//                System.out.println(coverageIncremental.peek().resultId+", "+getCurrentCost(resultPath)+", "+afterAddNewNodesPrice(resultPath, recordPaths, coverageIncremental.peek().resultId));
                while ((getCurrentCost(resultPath)+ afterAddNewNodesPrice(resultPath, recordPaths, coverageIncremental.peek().resultId))>budget){
//                    System.out.println("coverageIncremental.peek.resultId = "+coverageIncremental.peek().resultId + ",  "+ coverageIncremental.peek().lb+ ",  "+ coverageIncremental.peek().ub +",  "+afterAddNewNodesPrice(resultPath, recordPaths, coverageIncremental.peek().resultId));
                    coverageIncremental.poll();
//                    System.out.print(coverageIncremental.size()+", ");
                    if(coverageIncremental.isEmpty()){
                        break;
                    }
                } //find the maximum coverage incremental  and exceed the budget

                if (!coverageIncremental.isEmpty()){
                    relaxIndexNode currentOptimal = coverageIncremental.poll();
                    ArrayList<String> path = recordPaths.get(currentOptimal.resultId);
                    for (String tt: path){
//                        System.out.print(tt+", ");
                        resultPath.add(tt);
                    }

//                    if (maxRatio == false)
//                        System.out.println(currentOptimal.lb);
//                    else
//                        System.out.println(currentOptimal.ub);

                    for (Long l: S_L.get(currentOptimal.resultId)){
                        C.add(l);
                    }

                    leafNodeList.remove(currentOptimal.resultId);
//                    output
                    if (maxRatio == false){
                        rootleafNodeIncreHm.put(currentOptimal.resultId, (int) currentOptimal.lb);
//                        System.out.println("iteration "+iteration +" = "+currentOptimal.resultId +", coverage incremental = "+ currentOptimal.lb +", getCurrentCost(resultPath) = "+getCurrentCost(resultPath));
                    }else {
                        rootleafNodeIncreHm.put(currentOptimal.resultId, (int) currentOptimal.ub);
//                        System.out.println("iteration "+iteration +" = "+currentOptimal.resultId +", coverage incremental = "+ currentOptimal.ub +", getCurrentCost(resultPath) = "+getCurrentCost(resultPath));
                    }
                }else {
                    break;
                }
            }else {
                break;
            }
        }
//        System.out.println( "coverage = "+C.size()+ ",   getCurrentCost(resultPath) = "+ getCurrentCost(resultPath));
//        System.out.println();
        return C;

    }

    public static double assignPrice(){
//       先按照数据集占据的网格数分配价格，与网格数成正比,先假设价格就是  占据的网格数*0.1
//        datasetIdMapSignature.get(i).size()
        double totalPrice = 0.0;
        double unit = 0.10;
        double t = 0.80;
        for (String s: datasetIdMapSignature.keySet()){
//            t = t + unit;
//            String str22 = String.format("%.2f",t);
//            double price2 = Double.parseDouble(str22);
//            datasetPrice.put(s, price2);

            double price2 = datasetIdMapSignature.get(s).size()*unit;
            datasetPrice.put(s, price2); //price2

            totalPrice+= price2;
        }
//        datasetPrice.put("start", 0.0);
//        datasetPrice.put("end", 0.0);

        return totalPrice;
    }

    public static int findBest(HashMap<String, LinkedHashMap<String, Integer>> simpleGreedyResultS){
        int simpleMaxCoverage = 0;
        String simpleMaxId = "";
        for (String s1: simpleGreedyResultS.keySet()){
            int totalCoverage = 0;
            LinkedHashMap<String, Integer> lm = simpleGreedyResultS.get(s1);
            for (String s2: lm.keySet()){
                totalCoverage = totalCoverage+lm.get(s2);
            }
            if (simpleMaxCoverage<totalCoverage){
                simpleMaxCoverage = totalCoverage;
                simpleMaxId = s1;
            }
        }
        LinkedHashMap<String, Integer> simpleGreedyResult = simpleGreedyResultS.get(simpleMaxId);
//        output result
        for (String leafid: simpleGreedyResult.keySet()){
//            System.out.println(leafid+"," + integerMaptoDatasetId.get(Integer.parseInt(leafid))+",  coverage incremental = " +simpleGreedyResult.get(leafid)+", price = "+datasetPrice.get(leafid));
        }
        if (simpleMaxId.equals("")){
            System.out.print("simpleMaxId = -1"+",  "+simpleMaxCoverage+", ");

        }else {
            System.out.print("simpleMaxId = "+integerMaptoDatasetId.get(Integer.parseInt(simpleMaxId))+", "+simpleMaxCoverage+", ");

        }
        return simpleMaxCoverage;
    }
}

