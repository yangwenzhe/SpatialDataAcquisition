import Utils.indexAlgorithm;
import Utils.indexNode;

import java.io.*;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

public class DataLoader {
    public static double[] resolve(Long key, int resolution){
        double[] d = new double[2];
        d[0] = (int)(key % Math.pow(2,resolution)); //x
        d[1] = (int)(key / Math.pow(2,resolution));  //y
        return d;
    }

    public static long combine(double[] coordi, int resolution){
        long d = 0;
        double t = Math.pow(2,resolution);
        d = (long)(coordi[1]*t+coordi[0]);
//        for (int i=0; i<coordi.length-1; i++){
//            d += (coordi[i]*t);
//        }
//        d+= coordi[coordi.length-1];
        return d;
    }

    public static void writeFileAndMBR(HashMap<String, double[][]> zcodeMBR, String filePath, HashMap<String, HashMap<Long,Double>> dataSetMap) throws IOException {
        FileWriter fw = null;
        try {
            File file = new File(filePath);
            file.createNewFile();
            fw = new FileWriter(filePath);
            for (String id : dataSetMap.keySet()) {
                HashMap<Long,Double> map1 = dataSetMap.get(id);
                fw.write(id+",");
                double[][] MBRs = zcodeMBR.get(id);
                fw.write(MBRs[0][0]+",");
                fw.write(MBRs[0][1]+",");
                fw.write(MBRs[1][0]+",");
                fw.write(MBRs[1][1]+";");
//                fw.write(";");
                for (long key: map1.keySet()){
                    fw.write(String.valueOf(key));
                    fw.write(",");
                    String aa = String.valueOf(map1.get(key));
                    BigDecimal db = new BigDecimal(aa);
                    String ii = db.toPlainString();
                    fw.write(ii+";");
                }
                fw.write("\n");

            }
            fw.close();
            System.out.println("resolution.txt write successfully");

        }catch(Exception e){
            e.printStackTrace();
        }
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

            }
            fw.close();

        }catch(Exception e){
            e.printStackTrace();
        }
    }

    public static TreeMap<Long,Double> generateSignature(double minX, double minY, double rangeX, double rangeY, int resolution, String Path, String querySilo) throws CloneNotSupportedException{
        double blockSizeX = rangeX/Math.pow(2,resolution);
        double blockSizeY = rangeY/Math.pow(2,resolution);
//        TreeMap<Long,Integer> ZorderHis = new TreeMap<>();
        TreeMap<Long,Double> ZorderDensityHist = new TreeMap<>();
        long t = (long)Math.pow(2,resolution);
        int a,b;
        if (querySilo.equals("argoverse")){
            a = 6;
            b = 7;
        }else{
            a = 0;
            b = 1;
        }
        try {
            String record = "";
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(Path)));
            reader.readLine();
            while((record=reader.readLine())!=null){
                String[] fields = record.split(",");
                double x = Double.parseDouble(fields[a]) - minX;
                double y = Double.parseDouble(fields[b]) - minY;
                long X = (long) (x/blockSizeX);  //行row
                long Y = (long) (y/blockSizeY);  //列 col
                long id = Y*t+X; //Z-order code
                Double num = ZorderDensityHist.getOrDefault(id,new Double(0));
//                ZorderHis.put(id,num+1);
                ZorderDensityHist.put(id, (double)(num+1));
            }
            reader.close();
        }catch (IOException e) {
//            e.printStackTrace();
//            System.out.println("this dataset name isn't exist");
        }
        return ZorderDensityHist;
    }

    public static indexNode generateQuery(String dataPath, ArrayList<Long> arrayList, TreeMap<Long,Double> hashMap, int queryID, int dimension, int resolution, double minX,double minY,double rangeX, double rangeY, String querySilo) throws IOException,CloneNotSupportedException{
        //读查询数据集
        TreeMap<Long, Double> hashMap1= generateSignature(minX, minY, rangeX, rangeY, resolution, dataPath, querySilo);
        hashMap = (TreeMap<Long, Double>)hashMap1.clone();
        double[][] dataset = new double[hashMap.size()][dimension];
        int count = 0;
        System.out.print("queryID = "+queryID+".csv,  "+"size = "+hashMap.size()+"; ");
        for (long j: hashMap.keySet()){
            double[] d = DataLoader.resolve(j, resolution); // get coordinate
            dataset[count] = d;
            count++;
            arrayList.add(j);
        }
//        System.out.println();
        indexNode queryIndexNode = indexAlgorithm.createRootsDataset(dataset, dimension);
        queryIndexNode.rootToDataset = queryID;
        return queryIndexNode;
    }


    public static double distance2(double[] x, double[] y) {
        double d = 0.0;
        for (int i = 0; i < x.length; i++) {
            d += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return d;
    }
    public double distance(double[] x, double[] y) {
        return Math.sqrt(distance2(x, y));
    }


}
