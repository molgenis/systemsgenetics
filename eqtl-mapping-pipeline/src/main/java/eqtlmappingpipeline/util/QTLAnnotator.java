/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.regex.Pattern;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.probemapping.reading;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
class QTLAnnotator {

    public static void main(String[] args) throws IOException {
//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Great_GeneAssociations_20141127-public-2.0.2-aaQ50i-hg19-all-region_Ensembl.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLsFDR0.05-SNPLevel.txt",
//                "1;10-11;0;1;1;4-10", "probe;probe;probe", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered.txt-ExtendedInfoTesting.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved.txt",
//                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Normal\\RP3_2Mb_TSSextendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "1;1;1;10;1;11;1;4;1;4;1;4", "snp;probe;probe;probe;probe;probe", "D:\\UMCG\\ProbeMapping\\Info\\V70\\gencode.v15.annotation.gtf.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved-ExtendedInfo.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v5.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt", 
//                "1;8-9-10-11-12-13-14;1;4-5;0;17-18-20-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-70-71-81-133-134-141-142-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-192-193-194-195-196-197-198;0;1-2-3-4", "snp;probe;snp;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt-ExtendedInfo5.txt");
        
//            addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v10.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\statisticsTMM_exprssion.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSurounding_RP3Genes.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt",  
//                "0;8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40-41-42-43-44-45-46-47-48-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-69-70-71-72-73-74-75-76-77-78-79-80-81-82-83-84-85-86-87-88-89-90-91-92-93-94-95-96-97-98-99-100-101-102-103-104-105-106-107-108-109-110-111-112-113-114-115-116-117-118-119-120-121-122-123-124-125-126-127-128-129-130-131-132-133-134-135-136-137-138-139-140-141-142-143-144-145-146-147-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-192-193-194-195-196-197-198-199-200-201-202-203-204-205-206-207-208-209-210-211-212-213-214-215-216-217-218-219-220-221-222-223-224-225-226-227-228-229-230-231-232-233-234-235-236-237-238-239-240-241-242-243-244-245-246-247-248-249-250-251-252-253-254-255-256-257-258-259-260-261-262-263-264-265-266-267-268-269-270-271-272-273-274-275-276-277-278-279-280-281-282-283-284-285-286-287-288-289-290-291-292-293-294-295-296-297-298-299-300-301-302-303-304-305-306-307-308-309-310-311-312-313-314-315-316-317-318-319-320-321-322-323-324-325-326-327-328-329-330-331-332-333-334-335-336-337-338-339-340-341-342-343-344-345-346-347-348-349-350-351-352-353-354-355-356-357;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;1;4-5;0;1-2-3-4", "snp;snp;probe;probe;probe;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel_ExtendedInfTMP.txt");
//        
//            addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v10_13BM.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\statisticsTMM_exprssion.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSurounding_RP3Genes.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt",  
//                "0;8-9-10-11-12-13-14-15-16-17-18-175-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-192-193-194-195-196-197-198-199-200-201-202-203-204-205-206-207-208-291-292-293-294-295-296-297-298-299-300-301-302-303-304-305-306-307-308-309-310-311-312-313-314-315-316-317-318-319-320-321-322-323-324-325-326-327-328-329-330-331-332-333-334-335-336-337-338-339-340-341-342-343-344-345-346;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;1;4-5;0;1-2-3-4", "snp;snp;probe;probe;probe;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel_ExtendedInfo13BM.txt");
        
            addAnnotationToQTLOutput(
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_eQTMs_062015\\eQTLsFDR0.05-SNPLevel.txt",
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v10.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\statisticsTMM_exprssion.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSurounding_RP3Genes.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt",  
                "0;8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40-41-42-43-44-45-46-47-48-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-69-70-71-72-73-74-75-76-77-78-79-80-81-82-83-84-85-86-87-88-89-90-91-92-93-94-95-96-97-98-99-100-101-102-103-104-105-106-107-108-109-110-111-112-113-114-115-116-117-118-119-120-121-122-123-124-125-126-127-128-129-130-131-132-133-134-135-136-137-138-139-140-141-142-143-144-145-146-147-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-192-193-194-195-196-197-198-199-200-201-202-203-204-205-206-207-208-209-210-211-212-213-214-215-216-217-218-219-220-221-222-223-224-225-226-227-228-229-230-231-232-233-234-235-236-237-238-239-240-241-242-243-244-245-246-247-248-249-250-251-252-253-254-255-256-257-258-259-260-261-262-263-264-265-266-267-268-269-270-271-272-273-274-275-276-277-278-279-280-281-282-283-284-285-286-287-288-289-290-291-292-293-294-295-296-297-298-299-300-301-302-303-304-305-306-307-308-309-310-311-312-313-314-315-316-317-318-319-320-321-322-323-324-325-326-327-328-329-330-331-332-333-334-335-336-337-338-339-340-341-342-343-344-345-346-347-348-349-350-351-352-353-354-355-356-357;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;1;4-5;0;1-2-3-4", "snp;snp;probe;probe;probe;snp", null,
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_eQTMs_062015\\eQTLsFDR0.05-SNPLevel.txt_ExtendedInfo.txt");
        
//            addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Artificial_eQTMs0.0_Stringent.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v10.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\statisticsTMM_exprssion.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSurounding_RP3Genes.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt",  
//                "0;8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40-41-42-43-44-45-46-47-48-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-69-70-71-72-73-74-75-76-77-78-79-80-81-82-83-84-85-86-87-88-89-90-91-92-93-94-95-96-97-98-99-100-101-102-103-104-105-106-107-108-109-110-111-112-113-114-115-116-117-118-119-120-121-122-123-124-125-126-127-128-129-130-131-132-133-134-135-136-137-138-139-140-141-142-143-144-145-146-147-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-192-193-194-195-196-197-198-199-200-201-202-203-204-205-206-207-208-209-210-211-212-213-214-215-216-217-218-219-220-221-222-223-224-225-226-227-228-229-230-231-232-233-234-235-236-237-238-239-240-241-242-243-244-245-246-247-248-249-250-251-252-253-254-255-256-257-258-259-260-261-262-263-264-265-266-267-268-269-270-271-272-273-274-275-276-277-278-279-280-281-282-283-284-285-286-287-288-289-290-291-292-293-294-295-296-297-298-299-300-301-302-303-304-305-306-307-308-309-310-311-312-313-314-315-316-317-318-319-320-321-322-323-324-325-326-327-328-329-330-331-332-333-334-335-336-337-338-339-340-341-342-343-344-345-346-347-348-349-350-351-352-353-354-355-356-357;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15;1;4-5;0;1-2-3-4", "snp;snp;probe;probe;probe;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Artificial_eQTMs0.0_Stringent_ExtendedInfo.txt");

        
//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v5.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;",  
//                "1;8-9-10-11-12-13-14;1;4-5;0;17-18-20-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-70-71-81-133-134-141-142-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-192-193-194-195-196-197-198;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15", "snp;probe;snp;snp;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel_ExtendedInfo4.txt");
        
//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs_Exon\\Optimal_PC_and_QTL_Corrected\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ_v5.txt.gz;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\GC_ContentSuroundingProbes_full.txt;",  
//                "1;8-9-10-11-12-13-14;0;17-18-20-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-70-71-81-133-134-141-142-176-177-178-179-180-181-182-183-184-185-186-187-188-189-190-191-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174-175-192-193-194-195-196-197-198;0;1-2-3-4;0;1-2-3-4-5-6-7-8-9-10-11-12-13-14-15", "snp;snp;snp;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs_Exon\\Optimal_PC_and_QTL_Corrected\\eQTLSNPsFDR0.05-SNPLevel_ExtendedInfo2.txt");
//        
        
    }

    static void addAnnotationToQTLOutput(String in, String sources, String keyValuePairs, String idsToAnnotate, String reannotateGene, String out) throws IOException {
        QTLTextFile e = new QTLTextFile(in, QTLTextFile.R);

        ArrayList<EQTL> qtls = e.readList();
        String extraHeaders = "";
        ArrayList<String> extraAnnotation = new ArrayList<String>(qtls.size());

        String[] sourceList = sources.split(";");
        String[] keyValueList = keyValuePairs.split(";");
        String[] annotationIdList = idsToAnnotate.split(";");

        if (sourceList.length * 2 != keyValueList.length) {
            System.out.println("Faulty key value match input.");
            System.exit(0);
        }
        
        
        boolean addTo = false;
        
        if (!(reannotateGene == null || reannotateGene.equals(""))) {
            extraHeaders += "\tEnsemble Gene ids (;)";
            HashMap<String, String> gencodeGeneMapping = reading.readGTFAnnotationFileHash(reannotateGene, 50000);
            HashMap<String, HashSet<String>> genCodeMappingGenes = new HashMap<String, HashSet<String>>();
            for (Map.Entry<String, String> t : gencodeGeneMapping.entrySet()) {
                if (t.getKey().contains("ENSG")) {
                    String[] newKeys;
                    if (t.getValue().contains(";")) {
                        newKeys = t.getValue().split(";");
                    } else {
                        newKeys = new String[1];
                        newKeys[0] = t.getValue();
                    }

                    for (String nKey : newKeys) {
                        if (genCodeMappingGenes.containsKey(nKey)) {
                            genCodeMappingGenes.get(nKey).add(t.getKey());
                        } else {
                            HashSet<String> ensembleIds = new HashSet<String>();
                            ensembleIds.add(t.getKey());
                            genCodeMappingGenes.put(nKey, ensembleIds);
                        }
                    }
                }
            }
            int id = 0;
            for (EQTL qtl : qtls) {
                qtl.setProbeHUGO(qtl.getProbeHUGO().replace("\"", ""));
                
                if (qtl.getProbeHUGO().contains(";")) {
                    String[] ts = qtl.getProbeHUGO().split(";");
                    HashSet<String> setje = new HashSet<String>();
                    for (String x : ts) {
                        if (genCodeMappingGenes.containsKey(x)) {
                            setje.addAll(genCodeMappingGenes.get(x));
                        }
                    }
                    if(setje.isEmpty()){
                        if (!addTo) {
                            extraAnnotation.add("-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-" );
                        }
                    } else {
                        if (!addTo) {
                            extraAnnotation.add(setje.toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + setje.toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        }
                    }
                    
                } else {
                    if (genCodeMappingGenes.containsKey(qtl.getProbeHUGO())) {
                        if (!addTo) {
                            extraAnnotation.add(genCodeMappingGenes.get(qtl.getProbeHUGO()).toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + genCodeMappingGenes.get(qtl.getProbeHUGO()).toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        }
                    } else {
                        if (!addTo) {
                            extraAnnotation.add("-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-" );
                        }
                    }
                }
                id++;
            }
        }

        int masterId = 0;
        for (String annotationSource : sourceList) {
            if (!extraAnnotation.isEmpty()) {
                addTo = true;
            }
            if (annotationSource.contains("GWAS") && annotationSource.contains("Catalog")) {

                GWASCatalog cat = new GWASCatalog(annotationSource);
                HashMap<String, GWASSNP> catSnps = cat.getSnpToObj();
                extraHeaders += "\tGWASTraits (;)\tRisk allels (;)\t P-value GWAS (;)";

                int id = 0;
                for (EQTL eQtl : qtls) {
                    String gwasTrait = "";
                    String gwasAllelle = "";
                    String gwasPval = "";
                    GWASSNP s = catSnps.get(eQtl.getRsName());
                    if (s != null) {
                        HashSet<GWASTrait> traits = s.getAssociatedTraits();

                        for (GWASTrait t : traits) {
                            if (!gwasTrait.equals("")) {
                                gwasTrait += ";";
                                gwasAllelle += ";";
                                gwasPval += ";";
                            }
                            gwasTrait += t.getCleanName();
                            gwasAllelle += s.getRiskAllele(t);
                            gwasPval += s.getPValueAssociatedWithTrait(t);

                        }
//                        System.out.println(gwasTrait + "\t" + gwasAllelle+ "\t" + gwasPval);
                        if (!addTo) {
                            extraAnnotation.add(gwasTrait + "\t" + gwasAllelle + "\t" + gwasPval);
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + gwasTrait + "\t" + gwasAllelle + "\t" + gwasPval);
                        }


                    } else {
//                        System.out.println("-\t-\t-");
                        if (!addTo) {
                            extraAnnotation.add("-\t-\t-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-\t-\t-");
                        }

                    }
                    id++;
                }

            } else {
//                System.out.println(annotationSource);
//                System.out.println(masterId);
                int keyInt = 0;
                
                try {
                    keyInt = Integer.parseInt(keyValueList[(masterId * 2)]);
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }

                int[] valueInt = {0};
                try {
                    if(keyValueList[((masterId * 2) + 1)].contains("-")){
                        String[] tmpValues = keyValueList[((masterId * 2) + 1)].split("-");
                        valueInt = new int[tmpValues.length];
                        for(int i=0; i<valueInt.length; ++i){
                            valueInt[i] = Integer.parseInt(tmpValues[i]);
                        }
                    } 
                    else {
                        valueInt[0] = Integer.parseInt(keyValueList[((masterId * 2) + 1)]);
                    }
                    
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }
                
//                System.out.println(keyInt +"\t"+ valueInt);
                if(valueInt.length==1){
                    TextFile t = new TextFile(annotationSource, TextFile.R);
                    String[] header = t.readLineElems(Pattern.compile("\t"));
                    
                    extraHeaders += "\t" +header[valueInt[0]];

                    HashMap<String, String> annotationFromSource = (HashMap<String, String>) t.readAsHashMap(keyInt, valueInt[0]);

                    int id = 0;
                    for (EQTL eQtl : qtls) {

                        String newInfo;
                        if (annotationIdList[masterId].equalsIgnoreCase("probe")) {
                            newInfo = annotationFromSource.get(eQtl.getProbe());
                        } else {
                            newInfo = annotationFromSource.get(eQtl.getRsName());
                        }

                        if (newInfo != null) {
    //                        System.out.println(newInfo);
                            if (!addTo) {
                                extraAnnotation.add(newInfo);
                            } else {
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + newInfo);
                            }

                        } else {
                            System.out.println("Error, something wrong. Desired information is not pressent.");
                            System.out.println(annotationSource);
                            if (annotationIdList[masterId].equalsIgnoreCase("probe")) {
                                System.out.println(eQtl.getProbe());
                            } else {
                                System.out.println(eQtl.getRsName());
                            }
                            
//        System.exit(-1);
    //                        System.out.println("-");
                            if (!addTo) {
                                extraAnnotation.add("-");
                            } else {
                                
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t-");
                            }

                        }
                        id++;
                    }
                    if (!extraAnnotation.isEmpty()) {
                        addTo = true;
                    }
                    t.close();
                }  else if (valueInt.length>1) {
                    
                    TextFile t = new TextFile(annotationSource, TextFile.R);
                    String[] header = t.readLineElems(Pattern.compile("\t"));
                    for(int val : valueInt){
                        extraHeaders += "\t" +header[val];
                    }
                    

                    HashMap<String, String> annotationFromSource = (HashMap<String, String>) t.readAsHashMap(keyInt, valueInt, "\t");

                    int id = 0;
                    for (EQTL eQtl : qtls) {

                        String newInfo;
                        if (annotationIdList[masterId].equalsIgnoreCase("probe")) {
                            newInfo = annotationFromSource.get(eQtl.getProbe());
                        } else {
                            newInfo = annotationFromSource.get(eQtl.getRsName());
                        }

                        if (newInfo != null) {
    //                        System.out.println(newInfo);
                            if (!addTo) {
                                extraAnnotation.add(newInfo);
                            } else {
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + newInfo);
                            }

                        } else {
                            System.out.println("Error, something wrong. Desired information is not pressent.");
                            
                            if (annotationIdList[masterId].equalsIgnoreCase("probe")) {
                                System.out.println(eQtl.getProbe());
                            } else {
                                System.out.println(eQtl.getRsName());
                            }
                            
//        System.exit(-1);
//                            System.out.println("-");
                            if (!addTo) {
                                extraAnnotation.add("-");
                            } else {
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t-");
                            }

                        }
                        id++;
                    }
                    if (!extraAnnotation.isEmpty()) {
                        addTo = true;
                    }
                    t.close();
                }
            }
            masterId++;
        }

        TextFile outWriter = new TextFile(out, TextFile.W);

        outWriter.writeln(QTLTextFile.header + extraHeaders);
        int id = 0;
        for (EQTL qtl : qtls) {
            outWriter.writeln(qtl.toString() + "\t" + extraAnnotation.get(id));
            id++;
        }
        outWriter.close();
    }
}
