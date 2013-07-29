#Librerias importadas en la aplicacion
from Tkinter import*
from tkFileDialog import*
from tkColorChooser import*
from tkMessageBox import*
from rpy import*
from win32api import*
import math
import os
#Descripcion variables aplicacion
font_type=('Trebuchet MS',8)
font_type1=('Trebuchet MS',9)
font_type2=('Trebuchet MS',14,'bold','underline')
font_type3=('Trebuchet MS',11,'bold')
font_type4=('Trebuchet MS',10,'bold')
#Definicion dimensiones del Main Frame
win =Tk()
win.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
win.title('GeneQuery (Basic)')
width_resolution=[]
height_resolution=[]
win.resizable(0,0)
aff=[]
#nucleo principal del programa
def new_project():
    width_resolution.append(GetSystemMetrics (0))
    height_resolution.append(GetSystemMetrics (1))
    win.maxsize(width_resolution[len(width_resolution)-1],750)
    x_resolution=0
    y_resolution=0
    Frame(win,bg='#c9c9f1',width=1024,height=750,relief='groove').pack()
    #Listas definidas dentro de GeneQuery
    window_status=['close'];browse_status=['close'];up_window_status=['close'];user_workflow_status=['none'];help_status=['close'];data_transpose=['no']
    file_type=['none'];graph_type=['linear'];quality_list=['none'];annotate_list=[];active_workflow=['none'];unique_content=['BIOTIN']
    color=['#8080c0'];shape=['oval'];size=['15'];log_scale=['FALSE'];work_on_value=IntVar()
    design=[];samples=[];columns=[];control_samples=[];target_samples=[];folder=[];control_columns=[];file_name=[];microarray_list=['Affymetrix']
    values=[0];values2=[0];new_data2=[];norm_data1=[];new_data1=[];norm_data2=[];gene_id=[];data1=[];data2=[];tables=[]
    affymetrix_method=['RMA'];background_correction_method=['none'];analysis_type=['none'];genomic_technology=['none'];fdr_method_list=['close'];sam_statistic=['close'];baselines_list=['close']
    organisms_list=['close'];id_columns_list=['close'];id_types_list=['close'];similarities_list=['close'];clusters_list=['close'];chip_list=['close'];backgrounds_list=['close'];normalizations_list=['close'];imputate_list=['close'];aggregates_list=['close']
    column_list=['close'];row_list=['close'];value_list=['close'];aggregate_method=['none'];qpcr_missing=[35];missing_ct=StringVar()
    probes_list=['close'];column_normalization=IntVar();replaces=StringVar();analysis_method=IntVar();summarization_replaces=StringVar();replace_values=StringVar();fdr=StringVar();statistic=StringVar();graph_properties=StringVar()
    workflow=[];illumina_controls_file=[];position_x=[];position_y=[];scroll_position=[0]
    label_display=IntVar();graph_font=['Trebuchet MS'];graph_font_size=['6']
    graph_properties.set('Graph Options');analysis_method.set(1)
    pivot_columns=[];pivot_rows=[];pivot_datas=[]
    script_name=[];function_name=[];baseline_variable=[]
    intercept=[0];slope=[0];cor_value=[0];regression_window=['close'];tidy_up_mark=['FALSE']
    summarize=IntVar();impute=IntVar();normalize=IntVar();LIMMA=IntVar();SAM=IntVar();TTEST=IntVar();clustering=IntVar();biological_enrichment=IntVar()
    project=['close'];segmentation_list=[]
    col_name=['a','c','b','e','d','g','f','i','h','k','j','m','l','o','n','q','p','s','r','u','t','w','v','y','x','z','aa','cc','bb','ee','dd','gg','ff','ii','hh','kk','jj','mm','ll','oo','nn','qq','pp','ss','rr','uu','tt','ww','vv','yy','xx','zz']
    f=['a','c','b','e','d','g','f','i','h','k','j','m','l','o','n','q','p','s','r','u','t','w','v','y','x','z','aa','cc','bb','ee','dd','gg','ff','ii','hh','kk','jj','mm','ll','oo','nn','qq','pp','ss','rr','uu','tt','ww','vv','yy','xx','zz']
    background=StringVar();within_normalization=StringVar()
    #Interfaz grafica para el desarrollo de informes con los procesos analiticos usados
    def reporting():
        if analysis_type[len(analysis_type)-1]!='none':
            if window_status[len(window_status)-1]=='close':
                top=Toplevel(win)
                top.title ('GeneQuery report')
                top.resizable(0,0)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                background=['none'];normalize_correction=['none']
                Ct_change=['none'];aggregation=['none'];imputation=['none']
                technologies=['Affymetrix'];affymetrix_preprocess=StringVar();normalization_methods=StringVar();differential_expression_methods=StringVar();functional_analysis=IntVar();annotation_analysis=IntVar()
                def accept_report():
                    try:
                        file=open('C:\\PG\\processes\\data\\report.txt','w')
                        file.write('GeneQuery DATA ANALYSIS REPORT\n\n\n')
                        file.write('EXPERIMENT DETAILS\n\n')
                        file.write('Name:'+name_entry.get()+'\n')
                        file.write('Technology:'+technologies[len(technologies)-1]+'\n')
                        file.write('Analysis code:'+code_entry.get()+'\n\n')
                        file.write('Analysis summary:'+summary_entry.get(1.0,END)+'\n\n\n')
                        if technologies[len(technologies)-1]=='Affymetrix':
                            if affymetrix_preprocess.get()=='0':
                                preprocess='Data preprocess with MAS 5.0 method'
                            if affymetrix_preprocess.get()=='1':
                                preprocess='Data preprocess with RMA method'
                        if technologies[len(technologies)-1]=='Custom' or technologies[len(technologies)-1]=='Illumina':
                            preprocess='Data were not preprocessed with GeneQuery for genomics'
                        if technologies[len(technologies)-1]=='Agilent' or technologies[len(technologies)-1]=='Genepix' or technologies[len(technologies)-1]=='Quantarray' or technologies[len(technologies)-1]=='Scanarray' or technologies[len(technologies)-1]=='Imagene':
                            preprocess=background[len(background)-1]
                            within_normalization=normalize_correction[len(normalize_correction)-1]
                        file.write('ANALYTICAL FEATURES\n\n')
                        if technologies[len(technologies)-1]=='Affymetrix' or technologies[len(technologies)-1]=='Custom' or technologies[len(technologies)-1]=='Illumina':
                            file.write('Data preprocessing: '+preprocess+'\n\n')
                        if technologies[len(technologies)-1]=='Agilent' or technologies[len(technologies)-1]=='Genepix' or technologies[len(technologies)-1]=='Quantarray' or technologies[len(technologies)-1]=='Scanarray' or technologies[len(technologies)-1]=='Imagene':
                            file.write('Data preprocessing: '+background[len(background)-1]+' method was used for background correction. '+normalize_correction[len(normalize_correction)-1]+' method used for normalization within arrays\n\n')
                        if technologies[len(technologies)-1]=='qPCR':
                            file.write('Ct value for not detected genes:'+Ct_change[len(Ct_change)-1]+'\n')
                            file.write('Replicates aggregated by:'+aggregation[len(aggregation)-1]+'\n')
                            file.write('Imputation method applied:'+imputation[len(imputation)-1]+'\n\n')
                        if normalization_methods.get()=='0':
                            normalize_selected='Scale'
                        if normalization_methods.get()=='1':
                            normalize_selected='Quantile'
                        if normalization_methods.get()=='2':
                            normalize_selected='Z-Score'
                        if normalization_methods.get()=='3':
                            normalize_selected='Fold change'
                        file.write('Data Normalization: '+ 'Data were corrected by '+normalize_selected+' . Normalization consists on isolating statistical error in repeated measured data, this correction is necessary so that multiple chips can be compared to each other (making the distributions identical across arrays), and analyzed together.\n\n')
                        file.write('Differential expression analysis: The purpose of this step consists on study the differences in expression between two groups or classes of arrays, performing for each gene the corresponding test that will report a p-value.\n')
                        if differential_expression_methods.get()=='0':
                            file.write('Linear Models for Microarrays (LIMMA) was applied. LIMMA is a package for the analysis of gene expression microarray data, providing the ability to analyze comparisons between many RNA targets simultaneously in arbitrary complicated designed experiments. Empirical Bayesian methods are used to provide stable results even when the number of arrays is small.\n\n')
                        if differential_expression_methods.get()=='1':
                            file.write('Significant Analysis of Microarrays (SAM) was applied. SAM is a statistical technique for finding significant genes in a set of microarray experiments. The input to SAM consist on gene expression measurements from a set of microarray experiments, as well as a response variable from each experiment. The response variable may be a grouping like untreated, treated (either unpaired or paired), a multiclass grouping (like breast cancer, lymphoma, colon cancer, . . . ), a quantitative variable (like blood pressure) or a possibly censored survival time. SAM computes a statistic di for each gene i, measuring the strength of the relationship between gene expression and the response variable. It uses repeated permutations of the data to determine if the expression of any genes are significantly related to the response.\n\n')
                        if differential_expression_methods.get()=='2':
                            file.write('Unpaired T-Test was applied. Also known as Students t-test, it is applied to two different groups, assuming data are have got a normal distribution and standard deviation is approximately the same for both groups. The T-test compares the means of the two groups of data, determining whether the data come from the same population.\n\n')
                        if functional_analysis.get()==1:
                            file.write('Functional analysis: Aim of this process consists on the use of different sources of biological information to search for biological features (annotations) that frequently co-occur in a set of genes and rank them by statistical significance. It can be used to determine biological annotations or combinations of annotations that are significantly associated to a list of genes under study with respect to a reference list.\n\n')
                        if annotation_analysis.get()==1:
                            file.write('Annotation analysis: Descriptions, Gene sysmbols and other annotation IDs were added to the gene identifiers provided by using information contained within BioMart Database.')
                        file.close()
                        os.system('start winword c:\\PG\\processes\\data\\report.txt')
                        window_status.append('close')
                        top.destroy()
                    except:
                        warning=showinfo(title='Report',message='Impossible to connect with Microsoft Word')
                def quit_report():
                    top.destroy()
                    window_status.append('close')
                Frame(top, width=700, height=600,bg='#f5f5f5').pack()
                Frame(top, width=650,height=150,bg='#fcfcfc',borderwidth=2, relief='ridge').place(x=25,y=30)
                Frame(top, width=650,height=350,bg='#fcfcfc',borderwidth=2, relief='ridge').place(x=25,y=200)
                Label(top,text='Experiment details', font=font_type, bg='#f5f5f5', relief='groove').place(x=30,y=20)
                Label(top,text='Analytical features', font=font_type, bg='#f5f5f5', relief='groove').place(x=30,y=190)
                Label(top,text='Name:', font=font_type, bg='#fcfcfc').place(x=40,y=60)
                name_entry=Entry(top,width=30,font=font_type, bg='#fcfcfc',relief='groove')
                name_entry.place(x=120,y=60)
                def technology_selection():
                    def technology_selected(event):
                        index=technology_list.get(technology_list.curselection())
                        technology_list.destroy()
                        technology_list_scrollbar.destroy()
                        technology_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                        technology_entry.delete(0,END)
                        technology_entry.insert(INSERT,index)
                        technology_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                        technologies.append(index)
                        preprocessing_panel()
                        chip_list.append('close')
                    if chip_list[len(chip_list)-1]=='close':
                        technology_list=Listbox(top,width=30,height=6,font=font_type1,bg='#fcfcfc',relief='groove',cursor='hand2')
                        technology_list.place(x=120,y=90)
                        technology_list_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                        technology_list_scrollbar.place(x=305,y=90,height=120)
                        technology_list.config(yscrollcommand=technology_list_scrollbar.set)
                        technology_list_scrollbar.config(command=technology_list.yview)
                        technology_list.insert(END,'Affymetrix','Agilent','Custom','Genepix','Illumina','Imagene','qPCR','Quantarray','Scanarray')
                        technology_list.bind("<Double-Button-1>",technology_selected)
                        chip_list.append('open')
                Label(top,text='Technology:', font=font_type, bg='#fcfcfc').place(x=40,y=90)
                technology_entry=Entry(top,width=30,font=font_type, bg='#fcfcfc',relief='groove')
                technology_entry.place(x=120,y=90)
                technology_entry.insert(INSERT,'Affymetrix')
                technology_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',command=technology_selection,cursor='hand2').place(x=305,y=85)
                Label(top,text='Analysis code:', font=font_type, bg='#fcfcfc').place(x=40,y=120)
                code_entry=Entry(top,width=30,font=font_type, bg='#fcfcfc',relief='groove')
                code_entry.place(x=120,y=120)
                summary_entry=Text(top,width=50,height=4,font=font_type1, bg='#fcfcfc',relief='groove', wrap='word')
                summary_entry.place(x=350,y=65)
                Label(top,text='Summary', font=font_type, bg='#fcfcfc').place(x=355,y=50)
                def preprocessing_panel():
                    Frame(top,width=300, height=130,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=30,y=240)
                    Label(top,text='Data preprocessing', font=font_type1, bg='#fcfcfc').place(x=35,y=230)
                    if technologies[len(technologies)-1]=='Affymetrix':
                        Radiobutton(top,text='MAS 5.0',font=font_type,variable=affymetrix_preprocess,value=0,bg='#fcfcfc',cursor='hand2').place(x=35,y=260)
                        Radiobutton(top,text='RMA',font=font_type,variable=affymetrix_preprocess,value=1,bg='#fcfcfc',cursor='hand2').place(x=35,y=290)
                    if technologies[len(technologies)-1]=='Custom':
                        Label(top,text='Preprocessing not available for custom data',font=font_type,bg='#fcfcfc').place(x=35,y=260)
                    if technologies[len(technologies)-1]=='Illumina':
                        Label(top,text='Preprocessing not available for Illumina data',font=font_type,bg='#fcfcfc').place(x=35,y=260)
                    if technologies[len(technologies)-1]=='qPCR':
                        def aggregation_methods():
                            def aggregate_selection(event):
                                index=aggregation_list.get(aggregation_list.curselection())
                                aggregation_list.destroy()
                                aggregation_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                                aggregation_entry.delete(0,END)
                                aggregation_entry.insert(INSERT,index)
                                aggregation_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                                aggregates_list.append('close')
                            if aggregates_list[len(aggregates_list)-1]=='close':
                                aggregation_list=Listbox(top,width=15,height=2,font=font_type,bg='#fcfcfc', relief='groove')
                                aggregation_list.place(x=180,y=290)
                                aggregation_list.insert(END,'Median','Mean')
                                aggregation_list.bind("<Double-Button-1>",aggregate_selection)
                                aggregates_list.append('open')
                        def imputation_methods():
                            def imputate_selection(event):
                                index=imputation_list.get(imputation_list.curselection())
                                imputation_list.destroy()
                                imputation_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                                imputation_entry.delete(0,END)
                                imputation_entry.insert(INSERT,index)
                                imputation_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                                imputate_list.append('close')
                            if imputate_list[len(imputate_list)-1]=='close':
                                imputation_list=Listbox(top,width=15,height=3,font=font_type,bg='#fcfcfc', relief='groove')
                                imputation_list.place(x=180,y=320)
                                imputation_list.insert(END,'Constant','Row average','Column average')
                                imputation_list.bind("<Double-Button-1>",imputate_selection)
                                imputate_list.append('open')
                        Label(top,text='Ct for not detected genes:',font=font_type,bg='#fcfcfc').place(x=35,y=260)
                        change_Ct_entry=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
                        change_Ct_entry.place(x=180,y=260)
                        change_Ct_entry.insert(INSERT,'35')
                        Ct_change.append(change_Ct_entry.get())
                        Label(top,text='Replicates Aggregation:',font=font_type,bg='#fcfcfc').place(x=35,y=290)
                        aggregation_entry=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
                        aggregation_entry.place(x=180,y=290)
                        aggregation_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                        aggregation_entry.insert(INSERT,'Median')
                        aggregation.append(aggregation_entry.get())
                        aggregation_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                        Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=aggregation_methods).place(x=275,y=285)
                        Label(top,text='Imputation:',font=font_type,bg='#fcfcfc').place(x=35,y=320)
                        imputation_entry=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
                        imputation_entry.place(x=180,y=320)
                        imputation_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                        imputation_entry.insert(INSERT,'Constant')
                        imputation.append(imputation_entry.get())
                        imputation_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                        Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=imputation_methods).place(x=275,y=315)
                    if technologies[len(technologies)-1]=='Agilent' or technologies[len(technologies)-1]=='Genepix' or technologies[len(technologies)-1]=='Quantarray' or technologies[len(technologies)-1]=='Scanarray' or technologies[len(technologies)-1]=='Imagene':
                        def background_correction():
                            def background_selection(event):
                                index=background_list.get(background_list.curselection())
                                background_list.destroy()
                                background_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                                background_entry.delete(0,END)
                                background_entry.insert(INSERT,index)
                                background_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                                background.append(index)
                                backgrounds_list.append('close')
                            if backgrounds_list[len(backgrounds_list)-1]=='close':
                                background_list=Listbox(top,width=15,height=4,font=font_type,bg='#fcfcfc',relief='groove')
                                background_list.place(x=180,y=260)
                                background_list.insert(END,'none','subtract','minimum','normexp')
                                background_list.bind("<Double-Button-1>",background_selection)
                                backgrounds_list.append('open')
                        def within_normalization_methods():
                            def normalization_selection(event):
                                index=normalization_list.get(normalization_list.curselection())
                                normalization_list.destroy()
                                normalization_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                                normalization_entry.delete(0,END)
                                normalization_entry.insert(INSERT,index)
                                normalization_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                                normalize_correction.append(index)
                                normalizations_list.append('close')
                            if normalizations_list[len(normalizations_list)-1]=='close':
                                normalization_list=Listbox(top,width=15,height=3,font=font_type,bg='#fcfcfc',relief='groove')
                                normalization_list.place(x=180,y=290)
                                normalization_list.insert(END,'none','median','loess')
                                normalization_list.bind("<Double-Button-1>",normalization_selection)
                                normalizations_list.append('open')
                        Label(top,text='Background correction:',font=font_type,bg='#fcfcfc').place(x=35,y=260)
                        background_entry=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
                        background_entry.place(x=180,y=260)
                        background_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                        background_entry.insert(INSERT,'none')
                        background_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                        Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=background_correction).place(x=275,y=255)
                        Label(top,text='Within Normalization:',font=font_type,bg='#fcfcfc').place(x=35,y=290)
                        normalization_entry=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
                        normalization_entry.place(x=180,y=290)
                        normalization_entry.config(state=NORMAL,disabledbackground='#fcfcfc')
                        normalization_entry.insert(INSERT,'none')
                        normalization_entry.config(state=DISABLED,disabledbackground='#fcfcfc')
                        Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=within_normalization_methods).place(x=275,y=285)
                Frame(top,width=300, height=130,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=350,y=240)
                Label(top,text='Normalization', font=font_type, bg='#fcfcfc').place(x=355,y=230)
                Radiobutton(top,text='Scale',font=font_type,variable=normalization_methods,value=0,bg='#fcfcfc',cursor='hand2').place(x=360,y=260)
                Radiobutton(top,text='Quantile',font=font_type,variable=normalization_methods,value=1,bg='#fcfcfc',cursor='hand2').place(x=360,y=290)
                Radiobutton(top,text='Z-Score',font=font_type,variable=normalization_methods,value=2,bg='#fcfcfc',cursor='hand2').place(x=480,y=260)
                Radiobutton(top,text='Fold change',font=font_type,variable=normalization_methods,value=3,bg='#fcfcfc',cursor='hand2').place(x=480,y=290)
                Frame(top,width=300, height=130,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=30,y=390)
                Label(top,text='Differential expression', font=font_type, bg='#fcfcfc').place(x=35,y=380)
                Radiobutton(top,text='Linear Models for Microarrays',variable=differential_expression_methods,value=0,font=font_type,bg='#fcfcfc',cursor='hand2').place(x=40,y=410)
                Radiobutton(top,text='Significant Analysis of Microarrays',variable=differential_expression_methods,value=1,font=font_type,bg='#fcfcfc',cursor='hand2').place(x=40,y=440)
                Radiobutton(top,text='T-Test',variable=differential_expression_methods,value=2,font=font_type,bg='#fcfcfc',cursor='hand2').place(x=40,y=470)
                Frame(top,width=300, height=130,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=350,y=390)
                Label(top,text='Functional analysis', font=font_type, bg='#fcfcfc').place(x=355,y=380)
                Checkbutton(top,text='Gene annotation',font=font_type,variable=annotation_analysis,bg='#fcfcfc',cursor='hand2').place(x=360,y=410)
                Button(top,text='Help',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2').place(x=25,y=560)
                Button(top,text='OK',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2',command=accept_report).place(x=560,y=560)
                Button(top,text='Cancel',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2',command=quit_report).place(x=625,y=560)
                preprocessing_panel()
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',quit_report)
    #Definicion de las clases para realizar los analisis
    #ANALISIS FUNCIONAL
    class functional_analysis():
        def retrieve_annotations(self):
            if analysis_type.count('raw')>=1 or workflow.count('preprocessing')==1:
                if window_status[len(window_status)-1]=='close':
                    window_status.append('open')
                    organisms=['Homo sapiens'];organism_id=['hsapiens_gene_ensembl'];chip_id=['affy_hc_g110'];column_id=['GeneID'];annotations=[]
                    description=IntVar();embl=IntVar();ensembl=IntVar();entrezgene=IntVar();gene_id=IntVar();bp=IntVar();cc=IntVar();mf=IntVar()
                    hgnc=IntVar();interpro=IntVar();interpro_description=IntVar();merops=IntVar();mirbase=IntVar();pdb=IntVar();pfam=IntVar()
                    unigene=IntVar();uniprot=IntVar()
                    def cancel_database_annotation():
                        top.destroy()
                        window_status.append('close')
                        id_columns_list.append('close')
                        organisms_list.append('close')
                        id_types_list.append('close')
                    def accept_annotation():
                        del(annotate_list[0:len(annotate_list)])
                        if description.get()==1:
                            annotations.append('description')
                            annotate_list.append('true')
                        if embl.get()==1:
                            annotations.append('embl')
                            annotate_list.append('true')
                        if ensembl.get()==1:
                            annotations.append('ensembl_gene_id')
                            annotate_list.append('true')
                        if entrezgene.get()==1:
                            annotations.append('entrezgene')
                            annotate_list.append('true')
                        if gene_id.get()==1:
                            annotations.append('external_gene_id')
                            annotate_list.append('true')
                        if bp.get()==1:
                            annotations.append('go_biological_process_id')
                            annotate_list.append('true')
                        if cc.get()==1:
                            annotations.append('go_cellular_component_id')
                            annotate_list.append('true')
                        if mf.get()==1:
                            annotations.append('go_molecular_function_id')
                            annotate_list.append('true')
                        if hgnc.get()==1:
                            annotations.append('hgnc_id')
                            annotate_list.append('true')
                        if interpro.get()==1:
                            annotations.append('interpro')
                            annotate_list.append('true')
                        if interpro_description.get()==1:
                            annotations.append('interpro_description')
                            annotate_list.append('true')
                        if merops.get()==1:
                            annotations.append('merops')
                            annotate_list.append('true')
                        if mirbase.get()==1:
                            annotations.append('mirbase_id')
                            annotate_list.append('true')
                        if pdb.get()==1:
                            annotations.append('pdb')
                            annotate_list.append('true')
                        if pfam.get()==1:
                            annotations.append('pfam')
                            annotate_list.append('true')
                        if unigene.get()==1:
                            annotations.append('unigene')
                            annotate_list.append('true')
                        if uniprot.get()==1:
                            annotations.append('uniprot_swissprot')
                            annotate_list.append('true')
                        if len(annotations)>=1:
                            def fix_annotations():
                                up.destroy()
                                analysis_type.append('annotations')
                                message_label=Label(win,text='Fixing annotations...',font=font_type3,bg='#fcfcfc',relief='flat')
                                message_label.place(x=10,y=550)
                                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                                win.update()
                                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                                infile=r.matrix(infile)
                                ident=column_id[len(column_id)-1]
                                a=len(infile)
                                n=0
                                for i in range (0,a,1):
                                      if infile[i][0][0]==ident:
                                            n=i
                                identificador=[]
                                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                                win.update()
                                file=open('c:\\PG\\processes\\data\\result.txt','r')
                                content=file.readlines()
                                files=open('c:\\PG\\processes\\data\\annotation.txt','r')
                                infiles=files.readlines()
                                tab_number=(infiles[0].count('\t'))+1
                                for i in range(0,len(infiles),1):
                                      m=infiles[i].find('\t')
                                      ids=infiles[i][1:m-1]
                                      identificador.append(ids)
                                content[0]=content[0][0:len(content[0])-1]+'\t'+infiles[0]
                                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                                win.update()
                                duplicated_id=['ident']
                                for i in range (0,len(identificador),1):
                                      for j in range(0,len(infile[n][0]),1):
                                            if identificador[i]==infile[n][0][j]:
                                                    if duplicated_id.count(identificador[i])==0:
                                                        content[j]=content[j][0:len(content[j])-1]+'\t'+infiles[i]
                                                        duplicated_id.append(infile[n][0][j])
                                for i in range(0,len(content),1):
                                    if content[i].count('\t')!=content[0].count('\t'):
                                        content[i]=content[i][0:len(content[i])-1]+(tab_number*('\t'+'NA'))+'\n'
                                    content[i].replace("",'NA')
                                process.create_line(0,7,200,7,fill='#8080c0',width=15)
                                win.update()
                                result=open('c:\\PG\\processes\\data\\result.txt','w')
                                for i in range(0,len(content),1):
                                      result.write(content[i])
                                result.close()
                                process.delete(ALL)
                                message_label.destroy()
                                win.update()
                                data_table()
                            message_label=Label(win,text='Retrieving annotations...',font=font_type3,bg='#fcfcfc',relief='flat')
                            message_label.place(x=10,y=550)
                            process.create_line(0,7,50,7,fill='#8080c0',width=15)
                            win.update()
                            attributes=annotations[0:len(annotations)]
                            top.destroy()
                            r.source('c:\\PG\\srs.R')
                            process.create_line(0,7,100,7,fill='#8080c0',width=15)
                            win.update()
                            try:
                                r.databaseannotation('c:\\PG\\processes\\data\\result.txt',column_id[len(column_id)-1],organism_id[len(organism_id)-1],chip_id[len(chip_id)-1],attributes)
                                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                                win.update()
                                process.create_line(0,7,200,7,fill='#8080c0',width=15)
                                window_status.append('close')
                                win.update()
                                process.delete(ALL)
                                message_label.destroy()
                                infile=r.read_table('c:\\PG\\processes\\data\\annotation.txt',sep='\t')
                                infile=r.matrix(infile)
                                up=Toplevel(win)
                                up.title('Annotation results')
                                up.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                                up.maxsize(580,550)
                                up.resizable(0,0)
                                Frame(up,width=580,height=550).pack()
                                Frame(up,width=580,height=550,bg='#f5f5f5',relief='groove',borderwidth=2).place(x=0,y=0)
                                Label(up,text='Annotations associated to ID list',font=font_type3,bg='#f5f5f5').place(x=120,y=10)
                                rows=[]
                                for i in range(20):
                                    cols = []
                                    for j in range(5):
                                        e=Entry(up,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                                        e.place(x=10+j*109, y=45+20*i)
                                        cols.append(e)
                                    rows.append(cols)
                                b=infile
                                c=len(b[0][0])
                                d=len(b)
                                rows=[]
                                if c>20:
                                    c=20
                                if d>5:
                                    d=5
                                for i in range(c):
                                    cols = []
                                    for j in range(d):
                                        e=Entry(up,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                                        e.place(x=10+j*109, y=45+20*i)
                                        e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                        e.insert(INSERT,b[j][0][i])
                                        e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                        cols.append(e)
                                    rows.append(cols)
                                def close_annotation():
                                    up.destroy()
                                Button(up,text='Help',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=manuals).place(x=10,y=500)
                                Button(up,text='Fix table',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=fix_annotations).place(x=440,y=500)
                                Button(up,text='Close',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=close_annotation).place(x=510,y=500)
                                if active_workflow[len(active_workflow)-1]!='none':
                                    if genomic_technology[len(genomic_technology)-1]!='Custom':
                                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                            except:
                                warning=showinfo(title='Retrieve Annotations',message='Unable to get annotations requested')
                                message_label.destroy()
                                process.delete(ALL)
                                window_status.append('close')
                                if active_workflow[len(active_workflow)-1]!='none':
                                    if genomic_technology[len(genomic_technology)-1]!='Custom':
                                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)                                                                       
                        if len(annotations)==0:
                            warning=showinfo(title='Retrieve Annotations',message='One annotation term should be selected at least')
                    def organism_selection():
                            def organism_selected(event):
                                organisms_list.append('close')
                                index=organism_list.get(organism_list.curselection())
                                organism_list.destroy()
                                organism.config(state=NORMAL,disabledbackground='#ffffff')
                                organism.delete(0,END)
                                organism.insert(INSERT,index)
                                organism.config(state=DISABLED,disabledbackground='#ffffff')
                                organisms.append(index)
                                if index=='Arabidopsis thaliana':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_ath1_121501')
                                    organism_id.append('athaliana_eg_gene')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Homo sapiens':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_hc_g110')
                                    organism_id.append('hsapiens_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Mus musculus':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_mg_u74a')
                                    organism_id.append('mmusculus_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Rattus norvegicus':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_rae230a')
                                    organism_id.append('rnorvegicus_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Caenorhabditis elegans':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'agilent_oligoarray')
                                    organism_id.append('celegans_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Drosophila melanogaster':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_drosgenome1')
                                    organism_id.append('dmelanogaster_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                                if index=='Saccharomyces cerevisiae':
                                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                                    id_type.delete(0,END)
                                    id_type.insert(INSERT,'affy_yeast_2')
                                    organism_id.append('scerevisiae_gene_ensembl')
                                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                            if organisms_list[len(organisms_list)-1]=='close':
                                organisms_list.append('open')
                                organism_list=Listbox(top,width=40,height=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                                organism_list.place(x=80,y=20)
                                organism_list.insert(END,'Arabidopsis thaliana','Caenorhabditis elegans','Drosophila melanogaster','Homo sapiens','Mus musculus','Rattus norvegicus','Saccharomyces cerevisiae')
                                organism_list.bind("<Double-Button-1>",organism_selected)
                    def id_column_selection():
                        def id_column_selected(event):
                            index=id_column_list.get(id_column_list.curselection())
                            id_column_list.destroy()
                            id_column.config(state=NORMAL,disabledbackground='#ffffff')
                            id_column.delete(0,END)
                            id_column.insert(INSERT,index)
                            id_column.config(state=DISABLED,disabledbackground='#ffffff')
                            column_id.append(index)
                            id_column_scrollbar.destroy()
                            id_columns_list.append('close')
                        if id_columns_list[len(id_columns_list)-1]=='close':
                            id_columns_list.append('open')
                            id_column_list=Listbox(top,width=40,height=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                            id_column_list.place(x=80,y=50)
                            id_column_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                            id_column_scrollbar.place(x=325,y=50,height=123)
                            id_column_list.config(yscrollcommand=id_column_scrollbar.set)
                            id_column_scrollbar.config(command=id_column_list.yview)
                            if genomic_technology[len(genomic_technology)-1]=='Custom' and analysis_type[len(analysis_type)-1]=='raw':
                                for i in range(0,len(columns),1):
                                    id_column_list.insert(END,columns[i])
                            if genomic_technology[len(genomic_technology)-1]=='Custom' and analysis_type[len(analysis_type)-1]!='raw':
                                for i in range(0,len(columns),1):
                                    id_column_list.insert(END,columns[i])
                            if genomic_technology[len(genomic_technology)-1]!='Custom' and analysis_type[len(analysis_type)-1]!='raw':
                                for i in range(0,len(columns),1):
                                    id_column_list.insert(END,columns[i])
                            if genomic_technology[len(genomic_technology)-1]!='Custom' and analysis_type[len(analysis_type)-1]=='raw':
                                for i in range(0,len(columns)-1,1):
                                    id_column_list.insert(END,columns[i])
                            id_column_list.bind("<Double-Button-1>",id_column_selected)
                    def id_type_selection():
                        def id_type_selected(event):
                            index=id_type_list.get(id_type_list.curselection())
                            id_type_list.destroy()
                            id_type.config(state=NORMAL,disabledbackground='#ffffff')
                            id_type.delete(0,END)
                            id_type.insert(INSERT,index)
                            id_type.config(state=DISABLED,disabledbackground='#ffffff')
                            chip_id.append(index)
                            id_types_list.append('close')
                            id_type_scrollbar.destroy()
                        if id_types_list[len(id_types_list)-1]=='close':
                            id_types_list.append('open')
                            id_type_list=Listbox(top,width=40,height=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                            id_type_list.place(x=80,y=80)
                            id_type_list.bind("<Double-Button-1>",id_type_selected)
                            id_type_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                            id_type_scrollbar.place(x=325,y=80,height=123)
                            id_type_list.config(yscrollcommand=id_type_scrollbar.set)
                            id_type_scrollbar.config(command=id_type_list.yview)
                            if organisms[len(organisms)-1]=='Arabidopsis thaliana':
                                id_type_list.insert(END,'affy_ath1_121501','agilent_g2519f_015059','agilent_g2519f_021169','embl','ensembl_gene_id','entrezgene','go','interpro','uniprot_sptrembl')
                            if organisms[len(organisms)-1]=='Homo sapiens':
                                id_type_list.insert(END,'affy_hc_g110','affy_hg_focus','affy_hg_u133_plus_2','affy_hg_u133a','affy_hg_u133a_2','affy_hg_u133b','affy_hg_u95a','affy_hg_u95av2','affy_hg_u95b','affy_hg_u95c','affy_hg_u95d','affy_hg_u95e','affy_huex_1_0_st_v2','affy_hugene_1_0_st_v1','affy_u133_x3p','agilent_cgh_44b','efg_agilent_sureprint_g3_ge_8x60k','efg_agilent_wholegenome_4x44k_v1','efg_agilent_wholegenome_4x44k_v2','embl','ensembl_gene_id','entrezgene','go','hgcn_id','illumina_humanht_12','illumina_humanwg_6_v1','illumina_humanwg_6_v2','illumina_humanwg_6_v3','interpro','merops','mirbase_id','pdb','pfam','unigene','uniprot_swissprot')
                            if organisms[len(organisms)-1]=='Mus musculus':
                                id_type_list.insert(END,'affy_mg_u74a','affy_mg_u74av2','affy_mg_u74b','affy_mg_u74bv2','affy_mg_u74c','affy_mg_u74cv2','affy_moe430a','affy_moe430b','affy_moex_1_0_st_v1','affy_mogene_1_0_st_v1','affy_mouse430_2','affy_mouse430a_2','efg_agilent_sureprint_g3_ge_8x60k','efg_agilent_wholegenome_4x44k_v1','efg_agilent_wholegenome_4x44k_v2','embl','ensembl_gene_id','entrezgene','go','interpro','illumina_mousewg_6_v1','illumina_mousewg_6_v2','merops','mgi_id','mgi_symbol','mirbase_id','pdb','pfam','unigene','uniprot_swissprot')
                            if organisms[len(organisms)-1]=='Rattus norvegicus':
                                id_type_list.insert(END,'affy_rae230a','affy_rae230b','affy_raex_1_0_st_v1','affy_ragene_1_0_st_v1','affy_rat230_2','affy_rg_u34a','affy_rg_u34c','affy_rn_u34','affy_rt_u34','embl','entrezgene','go','interpro','pdb','pfam','rgd','rgd_symbol','unigene','uniprot_swissprot')
                            if organisms[len(organisms)-1]=='Caenorhabditis elegans':
                                id_type_list.insert(END,'agilent_oligoarray','embl','ensembl_gene_id','entrezgene','go','interpro','pdb','pfam','unigene','uniprot_swissprot')
                            if organisms[len(organisms)-1]=='Drosophila melanogaster':
                                id_type_list.insert(END,'affy_drosgenome1','affy_drosgenome2','embl','ensembl_gene_id','flybase_gene_id','flybasename_gene','go','interpro','mirbase_id','pdb','pfam','uniprot_swissprot')
                            if organisms[len(organisms)-1]=='Saccharomyces cerevisiae':
                                id_type_list.insert(END,'affy_yeast_2','affy_yg_s98','embl','ensembl_gene_id','entrezgene','go','interpro','sgd','uniprot_swissprot')
                    def check_all(event):
                        description.set(1);embl.set(1);ensembl.set(1);entrezgene.set(1);gene_id.set(1);bp.set(1);cc.set(1);mf.set(1);hgnc.set(1)
                        interpro.set(1);interpro_description.set(1);merops.set(1);mirbase.set(1);pdb.set(1);pfam.set(1);unigene.set(1);uniprot.set(1)
                    top=Toplevel(win)
                    top.title('Retrieve annotations')
                    top.resizable(0,0)
                    Frame(top, width=470,height=500,bg='#f5f5f5').pack()
                    top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                    Label(top,text='Organism:',font=font_type,bg='#f5f5f5',relief='flat').place(x=15,y=20)
                    Label(top,text='ID Column:',font=font_type,bg='#f5f5f5',relief='flat').place(x=15,y=50)
                    Label(top,text='ID Type:',font=font_type,bg='#f5f5f5',relief='flat').place(x=15,y=80)
                    organism=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
                    organism.place(x=80,y=20)
                    organism.config(state=NORMAL,disabledbackground='#ffffff')
                    organism.insert(INSERT,'Homo sapiens')
                    organism.config(state=DISABLED,disabledbackground='#ffffff')
                    id_column=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
                    id_column.place(x=80,y=50)
                    id_column.config(state=NORMAL,disabledbackground='#ffffff')
                    id_column.insert(INSERT,columns[0])
                    id_column.config(state=DISABLED,disabledbackground='#ffffff')
                    column_id.append(columns[0])
                    id_type=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
                    id_type.place(x=80,y=80)
                    id_type.config(state=NORMAL,disabledbackground='#ffffff')
                    id_type.insert(INSERT,'affy_hc_g110')
                    id_type.config(state=DISABLED,disabledbackground='#ffffff')
                    organism_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=organism_selection)
                    organism_button.place(x=325,y=15)
                    id_column_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=id_column_selection)
                    id_column_button.place(x=325,y=45)
                    id_type_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=id_type_selection)
                    id_type_button.place(x=325,y=75)
                    Frame(top,width=450,height=300,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=10,y=130)
                    Checkbutton(top,text='Description',font=font_type,variable=description,bg='#fcfcfc',cursor='hand2').place(x=20,y=150)
                    Checkbutton(top,text='EMBL',font=font_type,variable=embl,bg='#fcfcfc',cursor='hand2').place(x=20,y=180)
                    Checkbutton(top,text='ENSEMBL Gene ID',font=font_type,variable=ensembl,bg='#fcfcfc',cursor='hand2').place(x=20,y=210)
                    Checkbutton(top,text='Entrezgene',font=font_type,variable=entrezgene,bg='#fcfcfc',cursor='hand2').place(x=20,y=240)
                    Checkbutton(top,text='Gene symbol',font=font_type,variable=gene_id,bg='#fcfcfc',cursor='hand2').place(x=20,y=270)
                    Checkbutton(top,text='GO Biological Process',font=font_type,variable=bp,bg='#fcfcfc',cursor='hand2').place(x=20,y=300)
                    Checkbutton(top,text='GO Cellular Component',font=font_type,variable=cc,bg='#fcfcfc',cursor='hand2').place(x=20,y=330)
                    Checkbutton(top,text='GO Molecular Function ',font=font_type,variable=mf,bg='#fcfcfc',cursor='hand2').place(x=20,y=360)
                    Checkbutton(top,text='HGNC',font=font_type,variable=hgnc,bg='#fcfcfc',cursor='hand2').place(x=20,y=390)
                    Checkbutton(top,text='Interpro',font=font_type,variable=interpro,bg='#fcfcfc',cursor='hand2').place(x=220,y=150)
                    Checkbutton(top,text='Interpro Description',font=font_type,variable=interpro_description,bg='#fcfcfc',cursor='hand2').place(x=220,y=180)
                    Checkbutton(top,text='Merops',font=font_type,variable=merops,bg='#fcfcfc',cursor='hand2').place(x=220,y=210)
                    Checkbutton(top,text='Mirbase ID',font=font_type,variable=mirbase,bg='#fcfcfc',cursor='hand2').place(x=220,y=240)
                    Checkbutton(top,text='PDB',font=font_type,variable=pdb,bg='#fcfcfc',cursor='hand2').place(x=220,y=270)
                    Checkbutton(top,text='PFAM',font=font_type,variable=pfam,bg='#fcfcfc',cursor='hand2').place(x=220,y=300)
                    Checkbutton(top,text='Unigene',font=font_type,variable=unigene,bg='#fcfcfc',cursor='hand2').place(x=220,y=330)
                    Checkbutton(top,text='Uniprot',font=font_type,variable=uniprot,bg='#fcfcfc',cursor='hand2').place(x=220,y=360)
                    select_all=Label(top,text='Select all',font=font_type,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
                    select_all.place(x=240,y=390)
                    select_all.bind("<Button-1>",check_all)
                    Label(top,text='Annotations',font=font_type,bg='#f5f5f5',relief='ridge').place(x=15,y=120)
                    Button(top,text='Help',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2').place(x=10,y=450)
                    Button(top,text='Accept',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2',command=accept_annotation).place(x=340,y=450)
                    Button(top,text='Cancel',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2',command=cancel_database_annotation).place(x=405,y=450)
                    top.protocol('WM_DELETE_WINDOW',cancel_database_annotation)
    #CONTROLES DE CALIDAD
    class quality_class():
        def array_quality(self):
            def cancel_quality():
                window_status.append('close')
                top.destroy()
            if window_status[len(window_status)-1]=='close':
                top=Toplevel(win)
                top.resizable(0,0)
                top.title('Quality Controls')
                graph_colors=['#c9c9f1'];high_color=['#ff0000'];medium_color=['#000000'];low_color=['#00ff00'];graph_size=['9'];treatment=['raw']
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                window_status.append('open')
                control_label_size=StringVar()
                def graph_color_selection(event):
                    a=askcolor()
                    graph_colors.append(a[1])
                    graph_color=Button(top,width=7,bg=graph_colors[len(graph_colors)-1],relief='groove',cursor='hand2')
                    graph_color.place(x=180,y=360)
                    graph_color.bind("<Button-1>",graph_color_selection)
                def high_color_selection(event):
                    a=askcolor()
                    high_color.append(a[1])
                    high_label=Label(top,text='High',font=font_type,bg='#fcfcfc',fg=high_color[len(high_color)-1],cursor='hand2')
                    high_label.place(x=180,y=330)
                    high_label.bind("<Button-1>",high_color_selection)
                def medium_color_selection(event):
                    a=askcolor()
                    medium_color.append(a[1])
                    medium_label=Label(top,text='Medium',font=font_type,bg='#fcfcfc',fg=medium_color[len(medium_color)-1],cursor='hand2')
                    medium_label.place(x=220,y=330)
                    medium_label.bind("<Button-1>",medium_color_selection)
                def low_color_selection(event):
                    a=askcolor()
                    low_color.append(a[1])
                    low_label=Label(top,text='Low',font=font_type,bg='#fcfcfc',fg=low_color[len(low_color)-1],cursor='hand2')
                    low_label.place(x=280,y=330)
                    low_label.bind("<Button-1>",low_color_selection)
                def profile_view(event):
                    try:
                        if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                            r.source('C:\\PG\\srs.R')
                            r.profile('C:\\PG\\processes\\data\\hyb_result.txt')
                        if genomic_technology[len(genomic_technology)-1]!='Affymetrix' and genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                            warning=showinfo(title='Profiles',message='Visualization only available for Affymetrix data')
                    except:
                        pass
                def boxplot_view(event):
                    try:
                        font_size=int(control_label_size.get())
                        r.source('c:\\PG\\srs.R')
                        r.boxview('C:\\PG\\processes\\data\\raw_data.txt',graph_colors[len(graph_colors)-1],font_size)
                    except:
                        pass
                def normalized_boxplot_view(event):
                    try:
                        font_size=int(control_label_size.get())
                        if workflow.count('normalization')>0:
                            r.source('c:\\PG\\srs.R')
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                r.boxview('C:\\PG\\processes\\data\\Affymetrix_normalized_result.txt',graph_colors[len(graph_colors)-1],font_size)
                            if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Imagene' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                                r.boxview('C:\\PG\\processes\\data\\twochannel_normalized_result.txt',graph_colors[len(graph_colors)-1],font_size)
                            if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                                r.boxview('C:\\PG\\processes\\data\\normalized_result.txt',graph_colors[len(graph_colors)-1],font_size)
                        if workflow.count('normalization')==0:
                            warning=showinfo(title='Boxplot',message='Dataset not normalized yet')
                    except:
                        pass
                def correlation_map_view():
                    try:
                        r.source('c:\\PG\\srs.R')
                        if treatment[len(treatment)-1]=='raw':
                            r.circlecorr('C:\\PG\\processes\\data\\raw_data.txt',order='FALSE',col=(high_color[len(high_color)-1],medium_color[len(medium_color)-1],low_color[len(low_color)-1]),bg='gray50')
                        if treatment[len(treatment)-1]=='normalized':
                            if workflow.count('normalization')>0:
                                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                    r.circlecorr('C:\\PG\\processes\\data\\Affymetrix_normalized_result.txt',order='FALSE',col=(high_color[len(high_color)-1],medium_color[len(medium_color)-1],low_color[len(low_color)-1]),bg='gray50')
                                if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Imagene' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                                    r.circlecorr('C:\\PG\\processes\\data\\twochannel_normalized_result.txt',order='FALSE',col=(high_color[len(high_color)-1],medium_color[len(medium_color)-1],low_color[len(low_color)-1]),bg='gray50')
                                if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                                    r.circlecorr('C:\\PG\\processes\\data\\normalized_result.txt',order='FALSE',col=(high_color[len(high_color)-1],medium_color[len(medium_color)-1],low_color[len(low_color)-1]),bg='gray50')
                            if workflow.count('normalization')==0:
                                warning=showinfo(title='Boxplot',message='Dataset not normalized yet')
                    except:
                        pass
                def raw_correlation_view(event):
                    treatment.append('raw')
                    correlation_map_view()
                def normalized_correlation_view(event):
                    treatment.append('normalized')
                    correlation_map_view()
                def histogram_view(event):
                    def quit_histogram_view(event):
                        histogram_frame.destroy()
                        histogram_close.destroy()
                        sample_label.destroy()
                        sample_text.destroy()
                        scrollbar.destroy()
                    def sample_selection(event):
                        try:
                            index=sample_text.get(sample_text.curselection())
                            r.source('C:\\PG\\srs.R')
                            if genomic_technology[len(genomic_technology)-1]!='Illumina':
                                if genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                                    r.histogram('C:\\PG\\processes\\data\\result.txt',index,graph_colors[len(graph_colors)-1])
                                if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                                    r.histogram('C:\\PG\\processes\\data\\custom_result.txt',index,graph_colors[len(graph_colors)-1])    
                            if genomic_technology[len(genomic_technology)-1]=='Illumina':
                                r.histogram('C:\\PG\\processes\\data\\illumina_result.txt',index,graph_colors[len(graph_colors)-1])
                        except:
                            pass
                    histogram_frame=Frame(top,width=300,height=200,bg='#f5f5f5',borderwidth=2,relief='groove')
                    histogram_frame.place(x=50,y=100)
                    sample_text=Listbox(top,width=43,height=6,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                    sample_text.place(x=60,y=130)
                    sample_text.bind("<Double-Button-1>",sample_selection)
                    sample_label=Label(top,text='Sample Name',font=font_type,bg='#f5f5f5')
                    sample_label.place(x=60,y=105)
                    scrollbar=Scrollbar(top,bg='#f5f5f5')
                    scrollbar.place(x=325,y=130,height=115)
                    sample_text.config(yscrollcommand=scrollbar.set)
                    scrollbar.config(command=sample_text.yview)
                    for i in range (1,len(columns),1):
                        sample_text.insert(END,columns[i])
                    histogram_close=Button(top,text='Close',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2')
                    histogram_close.place(x=285,y=260)
                    histogram_close.bind("<Button-1>", quit_histogram_view)
                try:
                    Frame(top,width=400,height=420,bg='#fcfcfc').pack()
                    Label(top,image=foto29,bg='#fcfcfc').place(x=40,y=30)
                    Label(top,text='Control probes profiles',font=font_type,bg='#fcfcfc').place(x=80,y=30)
                    profile=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    profile.place(x=270,y=30)
                    profile.bind("<Button-1>",profile_view)
                    Label(top,image=foto31,bg='#fcfcfc').place(x=40,y=70)
                    Label(top,text='Correlation heatmap',font=font_type,bg='#fcfcfc').place(x=80,y=70)
                    correlation_view=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    correlation_view.place(x=270,y=70)
                    correlation_view.bind("<Button-1>", raw_correlation_view)
                    Label(top,image=foto13,bg='#fcfcfc').place(x=40,y=110)
                    Label(top,text='Raw data distribution',font=font_type,bg='#fcfcfc').place(x=80,y=110)
                    data_distribution=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    data_distribution.place(x=270,y=110)
                    data_distribution.bind("<Button-1>",histogram_view)
                    Label(top,image=foto32,bg='#fcfcfc').place(x=40,y=150)
                    Label(top,text='Raw data Boxplot',font=font_type,bg='#fcfcfc').place(x=80,y=150)
                    raw_boxplot_view=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    raw_boxplot_view.place(x=270,y=150)
                    raw_boxplot_view.bind("<Button-1>",boxplot_view)
                    Label(top,image=foto31,bg='#fcfcfc').place(x=40,y=190)
                    Label(top,text='Normalized correlation heatmap',font=font_type,bg='#fcfcfc').place(x=80,y=190)
                    normalized_correlation=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    normalized_correlation.place(x=270,y=190)
                    normalized_correlation.bind("<Button-1>",normalized_correlation_view)
                    Label(top,image=foto32,bg='#fcfcfc').place(x=40,y=230)
                    Label(top,text='Normalized data Boxplot',font=font_type,bg='#fcfcfc').place(x=80,y=230)
                    normalized_boxplot=Label(top,text='View',font=font_type,width=7,bg='#fcfcfc',fg='#0000ff',cursor='hand2')
                    normalized_boxplot.place(x=270,y=230)
                    normalized_boxplot.bind("<Button-1>",normalized_boxplot_view)
                    Frame(top,width=300,height=120,bg='#fcfcfc',borderwidth=2,relief='ridge').place(x=40,y=280)
                    Label(top,text='Settings',font=font_type,bg='#fcfcfc',relief='ridge',cursor='hand2').place(x=45,y=270)
                    Label(top,text='Label size:',font=font_type,bg='#fcfcfc').place(x=45,y=300)
                    label_size=Entry(top,width=2,font=font_type,bg='#fcfcfc',fg='#0000ff',textvariable=control_label_size,relief='flat',cursor='hand2')
                    label_size.place(x=180,y=300)
                    label_size.insert(INSERT,'9')
                    Label(top,text='Heatmap rank color:',font=font_type,bg='#fcfcfc').place(x=45,y=330)
                    high_label=Label(top,text='High',font=font_type,bg='#fcfcfc',fg=high_color[len(high_color)-1],cursor='hand2')
                    high_label.place(x=180,y=330)
                    high_label.bind("<Button-1>",high_color_selection)
                    medium_label=Label(top,text='Medium',font=font_type,bg='#fcfcfc',fg=medium_color[len(medium_color)-1],cursor='hand2')
                    medium_label.place(x=220,y=330)
                    medium_label.bind("<Button-1>",medium_color_selection)
                    low_label=Label(top,text='Low',font=font_type,bg='#fcfcfc',fg=low_color[len(low_color)-1],cursor='hand2')
                    low_label.place(x=280,y=330)
                    low_label.bind("<Button-1>",low_color_selection)
                    Label(top,text='Graph color:',font=font_type,bg='#fcfcfc').place(x=45,y=360)
                    graph_color=Button(top,width=7,bg=graph_colors[len(graph_colors)-1],relief='groove',cursor='hand2')
                    graph_color.place(x=180,y=360)
                    graph_color.bind("<Button-1>",graph_color_selection)
                    top.protocol('WM_DELETE_WINDOW',cancel_quality)
                    if genomic_technology[len(genomic_technology)-1]!='Custom':
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                except:
                    if genomic_technology[len(genomic_technology)-1]!='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=105)
    #PROCESAMIENTO Y CARGA DE DATOS 
    class preprocess_class():
        def data_pivoting(self):
            if genomic_technology[len(genomic_technology)-1]=='Custom':
                message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                message_label.place(x=10,y=550)
                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                win.update()
                r.source('C:\\PG\\srs.R')
                r.pivot('C:\\PG\\processes\\data\\pivot_result.txt',pivot_columns[len(pivot_columns)-1],pivot_rows[len(pivot_rows)-1],pivot_datas[len(pivot_datas)-1],'median')
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                file='C:\\PG\\processes\\data\\result.txt'
                content=open(file,'r')
                contents=content.readlines()
                column_number=(contents[0].count('\t'))+1
                infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
                infile=r.matrix(infile)
                del (columns[0:len(columns)])
                for i in range(0,len(infile),1):
                    columns.append(infile[i][0][0])
                for i in range(1,len(infile[1][0]),1):
                    gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                for i in range(0,len(data1),1):
                    if data1[i]>values[len(values)-1]:
                        values.append(data1[i])
                    if data2[i]>values2[len(values2)-1]:
                        values2.append(data2[i])
                    if values2[len(values2)-1]>values[len(values)-1]:
                        values.append(values2[len(values2)-1])
                for i in range(0,len(data1),1):
                    new_data1.append((data1[i]*650)/values[len(values)-1])
                    new_data2.append((data2[i]*400)/values[len(values)-1])
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                data_table()
                win.update()
                sample1.delete(0,END)
                sample2.delete(0,END)
                sample1.insert(INSERT,columns[2])
                sample2.insert(INSERT,columns[1])
                process.create_line(0,7,200,7,fill='#8080c0',width=15)
                win.update()
                graphic()
                process.delete(ALL)
                message_label.destroy()
        def hierarchical_clustering (self):
            if workflow.count('limma')>0 or analysis_type.count('sam')>0:
                del(samples[0:len(samples)])
                filter=IntVar();cluster_columns=IntVar();p_threshold=StringVar();adj_p_threshold=StringVar();cluster_method_selected=['complete'];similarity_method_selected=['euclidean']
                filter.set(2);cluster_high_color=['#ff0000'];cluster_medium_color=['#ffffff'];cluster_low_color=['#0000ff']
                def accept_hierarchical_clustering():
                    try:
                        workflow.append('cluster')
                        infile=open('c:\\PG\\processes\\data\\samples.txt','w')
                        for i in range(0,len(samples),1):
                            infile.write(samples[i]+'\n')
                        infile.close()
                        if len(samples)>1:
                            if samples.count('GeneID')==0 and samples.count('ID')==0 and samples.count('gene_id')==0:
                                if samples.count('logFC')>0 or samples.count('t')>0 or samples.count('P.Value')>0 or samples.count('adj.P.Val')>0 or samples.count('B')>0:
                                    warning=askquestion(title='',message='Some of the selected columns do not correspond to the Intensity signals.\nDo you want to continue?')
                                if samples.count('logFC')==0 and samples.count('t')==0 and samples.count('P.Value')==0 and samples.count('adj.P.Val')==0 and samples.count('B')==0:
                                    warning='yes'
                                if warning=='yes':                               
                                    analysis_type.append('cluster')
                                    top.destroy()
                                    window_status.append('close')
                                    win.update()
                                    if filter.get()==1:
                                        filter_type='P.Value'
                                        threshold_value=float(p_threshold.get())
                                    if filter.get()==2:
                                        filter_type='adj.P.Val'
                                        threshold_value=float(adj_p_threshold.get())
                                    if cluster_columns.get()==0:
                                        column='FALSE'
                                    if cluster_columns.get()==1:
                                        column='TRUE'
                                    process.create_line(0,7,50,7,fill='#8080c0',width=15)
                                    message_label=Label(win,text='CLUSTER analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                                    message_label.place(x=10,y=550)
                                    win.update()
                                    r.source('c:\\PG\\srs.R')
                                    process.create_line(0,7,100,7,fill='#8080c0',width=15)
                                    win.update()
                                    r.cluster('c:\\PG\\processes\\data\\result.txt',cluster_method_selected[len(cluster_method_selected)-1],similarity_method_selected[len(similarity_method_selected)-1],column,filter_type,threshold_value,cluster_high_color[len(cluster_high_color)-1],cluster_medium_color[len(cluster_medium_color)-1],cluster_low_color[len(cluster_low_color)-1])
                                    process.create_line(0,7,150,7,fill='#8080c0',width=15)
                                    win.update()
                                    if active_workflow[len(active_workflow)-1]!='none':
                                        if genomic_technology[len(genomic_technology)-1]!='Custom':
                                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=255)
                                        if genomic_technology[len(genomic_technology)-1]=='Custom':
                                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                                    information_text.config(state=NORMAL)
                                    information_text.insert(INSERT,'HIERARCHICAL CLUSTERING\n')
                                    information_text.insert(INSERT,'- Cluster method:'+cluster_method_selected[len(cluster_method_selected)-1]+'\n')
                                    information_text.insert(INSERT,'- Similarity method:'+similarity_method_selected[len(similarity_method_selected)-1]+'\n')
                                    information_text.insert(INSERT,'- Data filtered by '+filter_type+'\n')
                                    information_text.insert(INSERT,'- Threshold value: '+str(threshold_value)+'\n')
                                    information_text.config(state=DISABLED)
                                    win.update()
                                    message_label.destroy()
                                    process.delete(ALL)
                                    if warning=='no':
                                        pass
                            if samples.count('GeneID')>0 or samples.count('ID')>0 or samples.count('gene_id')>0:
                                warning=showinfo(title='Hierarchical Clustering',message='All columns selected must be numeric')
                        if len(samples)<=1:
                            warning=showinfo(title='Hierarchical Clustering',message='Two or more samples must be selected before Hierarchical Clustering analysis.')                            
                    except:
                        warning=showinfo(title='Hierarchical Clustering',message='Impossible to perform Hierarchical clustering. Review settings defined.')
                        message_label.destroy()
                        process.delete(ALL)
                        if active_workflow[len(active_workflow)-1]!='none':
                            if genomic_technology[len(genomic_technology)-1]!='Custom':
                                Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=255)
                            if genomic_technology[len(genomic_technology)-1]=='Custom':
                                Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                def cluster_selection(event):
                    def cluster_method_selection(event):
                        index=cluster_list.get(cluster_list.curselection())
                        cluster_method_selected.append(index)
                        cluster_method.config(state=NORMAL,disabledbackground='#ffffff')
                        cluster_method.delete(0,END)
                        cluster_method.insert(INSERT,index)
                        cluster_method.config(state=DISABLED,disabledbackground='#ffffff')
                        cluster_list.destroy()
                        clusters_list.append('close')
                    if clusters_list[len(clusters_list)-1]=='close':
                        clusters_list.append('open')
                        cluster_list=Listbox(top,width=20,height=5,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                        cluster_list.place(x=120,y=250)
                        cluster_list.insert(END,'complete','single', 'average', 'median', 'centroid')
                        cluster_list.bind("<Double-Button-1>",cluster_method_selection)
                def similarity_selection(event):
                    def similarity_method_selection(event):
                        index=similarity_list.get(similarity_list.curselection())
                        similarity_method_selected.append(index)
                        similarity_measure.config(state=NORMAL,disabledbackground='#ffffff')
                        similarity_measure.delete(0,END)
                        similarity_measure.insert(INSERT, index)
                        similarity_measure.config(state=DISABLED,disabledbackground='#ffffff')
                        similarity_list.destroy()
                        similarities_list.append('close')
                    if similarities_list[len(similarities_list)-1]=='close':
                        similarities_list.append('open')
                        similarity_list=Listbox(top,width=20,height=4,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                        similarity_list.place(x=120,y=290)
                        similarity_list.insert(END,'euclidean','maximum','manhattan','pearson')
                        similarity_list.bind("<Double-Button-1>",similarity_method_selection)
                def select_columns(event):
                    try:
                        indexes=available_columns.curselection()
                        for i in range (0,len(indexes),1):
                            index=available_columns.get(indexes[i])
                            selected_columns.insert(END,index)
                            samples.append(index)
                        for i in range (0,len(indexes),1):
                            available_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                def deselect_columns(event):
                    try:
                        indexes=selected_columns.curselection()
                        for i in range(0,len(indexes),1):
                            index=selected_columns.get(indexes[i])
                            if len(samples)==1:
                                del(samples[0])
                            for i in range (0,len(samples)-1,1):
                                if samples[i]==index:
                                    del(samples[i])
                            available_columns.insert(END,index)
                        for i in range(0,len(indexes),1):
                            selected_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                if window_status[len(window_status)-1]=='close':                    
                    def close_window():
                        top.destroy()
                        window_status.append('close')
                    def cluster_high_color_selection(event):
                        a=askcolor()
                        cluster_high_color.append(a[1])
                        high_label=Label(top,text='High',font=font_type,bg='#f5f5f5',fg=cluster_high_color[len(cluster_high_color)-1],cursor='hand2')
                        high_label.place(x=130,y=380)
                        high_label.bind("<Button-1>",cluster_high_color_selection)
                    def cluster_medium_color_selection(event):
                        a=askcolor()
                        cluster_medium_color.append(a[1])
                        medium_label=Label(top,text='Medium',font=font_type,bg='#f5f5f5',fg=cluster_medium_color[len(cluster_medium_color)-1],cursor='hand2')
                        medium_label.place(x=180,y=380)
                        medium_label.bind("<Button-1>",cluster_medium_color_selection)
                    def cluster_low_color_selection(event):
                        a=askcolor()
                        cluster_low_color.append(a[1])
                        low_label=Label(top,text='Low',font=font_type,bg='#f5f5f5',fg=cluster_low_color[len(cluster_low_color)-1],cursor='hand2')
                        low_label.place(x=250,y=380)
                        low_label.bind("<Button-1>",cluster_low_color_selection)
                    window_status.append('open')
                    top=Toplevel(win)
                    top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                    top.title('Hierarchical Clustering')
                    top.resizable(0,0)
                    top.maxsize(520,460)
                    Frame(top,bg='#f5f5f5',width=520,height=460,relief='groove').pack()
                    Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                    file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                    file_entry.place(x=80,y=10)
                    Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                    Label(top,text='Selected Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=300,y=40)
                    available_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                    available_columns.place(x=10,y=60)
                    available_columns.bind("<Double-Button-1>",select_columns)
                    available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    available_columns_scrollbar.place(x=190,y=60,height=175)
                    available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                    available_columns_scrollbar.config(command=available_columns.yview)
                    selected_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                    selected_columns.place(x=300,y=60)
                    selected_columns.bind("<Double-Button-1>",deselect_columns)
                    selected_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    selected_columns_scrollbar.place(x=480,y=60,height=175)
                    selected_columns.config(yscrollcommand=selected_columns_scrollbar.set)
                    selected_columns_scrollbar.config(command=selected_columns.yview)
                    accept=Button(top,image=foto26,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    accept.place(x=230,y=110)
                    accept.bind("<Button-1>",select_columns)
                    accept.bind("<Double-Button-1>",select_columns)
                    remove=Button(top,image=foto27,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    remove.place(x=230,y=150)
                    remove.bind("<Button-1>",deselect_columns)
                    remove.bind("<Double-Button-1>",deselect_columns)
                    Label(top,text='Clustering method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=250)
                    Label(top,text='Similarity measure:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=290)
                    cluster_method=Entry(top,width=20,bg='#ffffff',font=font_type,relief='groove')
                    cluster_method.place(x=120,y=250)
                    cluster_method.config(state=NORMAL,disabledbackground='#ffffff')
                    cluster_method.insert(INSERT,'complete')
                    cluster_method.config(state=DISABLED,disabledbackground='#ffffff')
                    cluster_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    cluster_button.place(x=245,y=245)
                    cluster_button.bind("<Button-1>",cluster_selection)
                    similarity_measure=Entry(top,width=20,bg='#ffffff',font=font_type,relief='groove')
                    similarity_measure.place(x=120,y=290)
                    similarity_measure.config(state=NORMAL,disabledbackground='#ffffff')
                    similarity_measure.insert(INSERT,'euclidean')
                    similarity_measure.config(state=DISABLED,disabledbackground='#ffffff')
                    similarity_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    similarity_button.place(x=245,y=285)
                    similarity_button.bind("<Button-1>",similarity_selection)
                    Frame(top,width=185,height=140,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=300,y=260)
                    Checkbutton(top,text='Reorder columns?',font=font_type,bg='#f5f5f5',variable=cluster_columns,relief='flat',cursor='hand2').place(x=10,y=340)
                    Label(top,text='Heatmap rank color:',font=font_type,bg='#f5f5f5').place(x=10,y=380)
                    cluster_high=Label(top,text='High',font=font_type,bg='#f5f5f5',fg=cluster_high_color[len(cluster_high_color)-1],cursor='hand2')
                    cluster_high.place(x=130,y=380)
                    cluster_high.bind("<Button-1>",cluster_high_color_selection)
                    cluster_medium=Label(top,text='Medium',font=font_type,bg='#f5f5f5',fg=cluster_medium_color[len(cluster_medium_color)-1],cursor='hand2')
                    cluster_medium.place(x=180,y=380)
                    cluster_medium.bind("<Button-1>",cluster_medium_color_selection)
                    cluster_low=Label(top,text='Low',font=font_type,bg='#f5f5f5',fg=cluster_low_color[len(cluster_low_color)-1],cursor='hand2')
                    cluster_low.place(x=250,y=380)
                    cluster_low.bind("<Button-1>",cluster_low_color_selection)
                    Label(top,text='Work on',font=font_type,bg='#f5f5f5',relief='ridge').place(x=310,y=250)
                    Radiobutton(top,text='P.Value<=',font=font_type,bg='#fcfcfc',variable=filter,value=1,relief='flat',cursor='hand2').place(x=310,y=290)
                    p_value=Entry(top,width=7,font=font_type,bg='#fcfcfc',textvariable=p_threshold,relief='groove')
                    p_value.place(x=410,y=290)
                    p_value.insert(INSERT,'0.05')
                    Radiobutton(top,text='Adj. P.Value<=',font=font_type,bg='#fcfcfc',variable=filter,value=2,relief='flat',cursor='hand2').place(x=310,y=320)
                    adj_p_value=Entry(top,width=7,font=font_type,bg='#fcfcfc',textvariable=adj_p_threshold,relief='groove')
                    adj_p_value.place(x=410,y=320)
                    adj_p_value.insert(INSERT,'0.05')
                    Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_hierarchical_clustering).place(x=370,y=420)
                    Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=close_window).place(x=425,y=420)
                    Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                    top.protocol('WM_DELETE_WINDOW',close_window)
                    for i in range(0,len(columns),1):
                        available_columns.insert(END,columns[i])
            if workflow.count('limma')==0 and analysis_type.count('sam')==0:
                warning=showinfo(title='GeneQuery Cluster',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
        def custom_script(self):
            if len(file_name)>0:
                if window_status[len(window_status)-1]=='close':
                    def close_window():
                        top.destroy()
                        window_status.append('close')
                    def select_columns(event):
                        try:
                            indexes=available_columns.curselection()
                            for i in range (0,len(indexes),1):
                                index=available_columns.get(indexes[i])
                                selected_columns.insert(END,index)
                                samples.append(index)
                            for i in range (0,len(indexes),1):
                                available_columns.delete(int(indexes[i])-i)
                        except:
                            warning=showwarning(title='Select Columns',message='No column selected')
                    def deselect_columns(event):
                        try:
                            indexes=selected_columns.curselection()
                            for i in range(0,len(indexes),1):
                                index=selected_columns.get(indexes[i])
                                if len(samples)==1:
                                    del(samples[0])
                                for i in range (0,len(samples)-1,1):
                                    if samples[i]==index:
                                        del(samples[i])
                                available_columns.insert(END,index)
                            for i in range(0,len(indexes),1):
                                selected_columns.delete(int(indexes[i])-i)
                        except:
                            warning=showwarning(title='Select Columns',message='No column selected')
                    def script_label(event):
                        def quit_script_label(event):
                            show_label.destroy()
                        show_label=Label(top,text='Browse Script',font=font_type,bg='#fcfcfc',relief='groove')
                        show_label.place(x=305,y=305)
                        script_button.bind("<Leave>",quit_script_label)
                    def browse_script(event):
                        file=askopenfilename()
                        script_entry.delete(0,END)
                        script_entry.insert(INSERT,file)
                    def accept_custom_script():
                        if len(samples)>1:
                            if samples.count('GeneID')==0 and samples.count('ID')==0 and samples.count('gene_id')==0:
                                script_name.append(script_entry.get())
                                function_name.append(function_entry.get())
                                workflow.append('script')
                                analysis_type.append('script')
                                top.destroy()
                                window_status.append('close')
                                data_table()
                            if samples.count('GeneID')>0 or samples.count('ID')>0 or samples.count('gene_id')>0:
                                warning=showinfo(title='Custom Script',message='All columns selected must be numeric')
                                window_status.append('close')
                                top.destroy()
                        if len(samples)<=1:
                            warning=showinfo(title='Custom Script',message='Two or more samples must be selected before using selected script.')
                    top=Toplevel(win)
                    top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                    top.title('Custom Script')
                    top.resizable(0,0)
                    top.maxsize(520,460)
                    Frame(top,bg='#f5f5f5',width=520,height=460,relief='groove').pack()
                    Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                    file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                    file_entry.place(x=80,y=10)
                    file_entry.config(state=NORMAL,disabledbackground='#f5f5f5')
                    if analysis_type[len(analysis_type)-1]=='Affymetrix' or analysis_type[len(analysis_type)-1]=='All-Affymetrix':
                        file_entry.insert(INSERT,'c:\\PG\\processes\\data\\result.txt')
                    if analysis_type[len(analysis_type)-1]!='Affymetrix':
                        file_entry.insert(INSERT,file_name[len(file_name)-1])
                    file_entry.config(state=DISABLED,disabledbackground='#f5f5f5')
                    Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                    Label(top,text='Selected Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=300,y=40)
                    available_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                    available_columns.place(x=10,y=60)
                    available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    available_columns_scrollbar.place(x=190,y=60,height=175)
                    available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                    available_columns_scrollbar.config(command=available_columns.yview)
                    selected_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                    selected_columns.place(x=300,y=60)
                    selected_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    selected_columns_scrollbar.place(x=480,y=60,height=175)
                    selected_columns.config(yscrollcommand=selected_columns_scrollbar.set)
                    selected_columns_scrollbar.config(command=selected_columns.yview)
                    available_columns.bind("<Double-Button-1>",select_columns)
                    selected_columns.bind("<Double-Button-1>",deselect_columns)
                    accept=Button(top,image=foto26,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    accept.place(x=230,y=110)
                    accept.bind("<Button-1>",select_columns)
                    accept.bind("<Double-Button-1>",select_columns)
                    remove=Button(top,image=foto27,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    remove.place(x=230,y=150)
                    remove.bind("<Button-1>",deselect_columns)
                    remove.bind("<Double-Button-1>",deselect_columns)
                    Label(top,text='Script file:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=280)
                    script_entry=Entry(top,width=30,font=font_type,bg='#ffffff',relief='groove')
                    script_entry.place(x=100,y=280)
                    script_button=Button(top,image=foto1,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                    script_button.place(x=290,y=275)
                    script_button.bind("<Enter>",script_label)
                    script_button.bind("<Button-1>",browse_script)
                    Label(top,text='Function name:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=320)
                    function_entry=Entry(top,width=30,font=font_type,bg='#ffffff',relief='groove')
                    function_entry.place(x=100,y=320)
                    Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_custom_script).place(x=370,y=420)
                    Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=close_window).place(x=425,y=420)
                    Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                    window_status.append('open')
                    for i in range(0,len(columns),1):
                        available_columns.insert(END,columns[i])
                    top.protocol('WM_DELETE_WINDOW',close_window)
            if len(file_name)==0:
                warning=showinfo(title='Custom Script',message='Dataset must be imported before using the Custom Script tool')
        def data_normalization(self):
            del(samples[0:len(samples)])
            def normalization_cancel():
                top.destroy()
                window_status.append('close')
                win.update()
            def select_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        selected_columns.insert(END,index)
                        samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_columns(event):
                try:
                    indexes=selected_columns.curselection()
                    for i in range(0,len(indexes),1):
                        index=selected_columns.get(indexes[i])
                        if len(samples)==1:
                            del(samples[0])
                        for i in range (0,len(samples)-1,1):
                            if samples[i]==index:
                                del(samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        selected_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def baseline_samples(event):
                def select_baseline(event):
                    index=baseline_list.get(baseline_list.curselection())
                    baseline.delete(0,END)
                    baseline.insert(INSERT,index)
                    baseline_list.destroy()
                    baselines_list.append('close')
                if baselines_list[len(baselines_list)-1]=='close':
                    baselines_list.append('open')
                    baseline_list=Listbox(top,width=9,height=3,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                    baseline_list.place(x=390,y=330)
                    for i in range(1,len(columns),1):
                        baseline_list.insert(END,columns[i])
                    baseline_list.bind("<Double-Button-1>",select_baseline)
            def accept_normalization():
                if len(samples)>1:
                    if samples.count('GeneID')==0 and samples.count('ID')==0 and samples.count('gene_id')==0:
                        workflow.append('normalization')
                        analysis_type.append('normalization')
                        baseline_variable.append(baseline.get())
                        top.destroy()
                        window_status.append('close')
                        top.update()
                        data_table()
                    if samples.count('GeneID')>0 or samples.count('ID')>0 or samples.count('gene_id')>0:
                        warning=showinfo(title='Normalization',message='All columns selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if len(samples)<=1:
                    warning=showinfo(title='Normalization',message='Two or more samples must be selected before data normalization.')
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Normalization')
                top.resizable(0,0)
                top.maxsize(520,460)
                Frame(top,bg='#f5f5f5',width=520,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                if analysis_type[len(analysis_type)-1]=='Affymetrix' or analysis_type[len(analysis_type)-1]=='All-Affymetrix':
                    file_entry.insert(INSERT,'c:\\PG\\processes\\data\\result.txt')
                if analysis_type[len(analysis_type)-1]!='Affymetrix' or analysis_type[len(analysis_type)-1]!='All-Affymetrix':
                    file_entry.insert(INSERT,file_name[len(file_name)-1])
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                Label(top,text='Selected Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=300,y=40)
                available_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=190,y=60,height=175)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                selected_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2', relief='groove',selectmode=MULTIPLE)
                selected_columns.place(x=300,y=60)
                selected_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                selected_columns_scrollbar.place(x=480,y=60,height=175)
                selected_columns.config(yscrollcommand=selected_columns_scrollbar.set)
                selected_columns_scrollbar.config(command=selected_columns.yview)
                available_columns.bind("<Double-Button-1>",select_columns)
                selected_columns.bind("<Double-Button-1>",deselect_columns)
                Frame(top,bg='#fcfcfc',width=185,height=150,relief='groove',borderwidth=2).place(x=10,y=250)
                Frame(top,bg='#fcfcfc',width=185,height=150,relief='groove',borderwidth=2).place(x=300,y=250)
                Label(top,text='Column Normalization',font=font_type,bg='#f5f5f5',relief='ridge').place(x=13,y=240)
                Label(top,text='Row Normalization',font=font_type,bg='#f5f5f5',relief='ridge').place(x=303,y=240)
                accept=Button(top,image=foto26,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                accept.place(x=230,y=110)
                accept.bind("<Button-1>",select_columns)
                accept.bind("<Double-Button-1>",select_columns)
                remove=Button(top,image=foto27,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                remove.place(x=230,y=150)
                remove.bind("<Button-1>",deselect_columns)
                remove.bind("<Double-Button-1>",deselect_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_normalization).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=normalization_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                Radiobutton(top,text='Scale',font=font_type,bg='#fcfcfc',variable=column_normalization,value=0,cursor='hand2').place(x=20,y=270)
                Radiobutton(top,text='Normalize by quantiles',font=font_type,bg='#fcfcfc',variable=column_normalization,value=1,cursor='hand2').place(x=20,y=300)
                Radiobutton(top,text='Z-Score',font=font_type,bg='#fcfcfc',variable=column_normalization,value=2,cursor='hand2').place(x=310,y=270)
                Radiobutton(top,text='Fold change',font=font_type,bg='#fcfcfc',variable=column_normalization,value=3,cursor='hand2').place(x=310,y=300)
                Label(top,text='Baseline value:',font=font_type,bg='#fcfcfc',relief='flat').place(x=310,y=330)
                baseline=Entry(top,width=9,font=font_type,bg='#fcfcfc',relief='groove')
                baseline.place(x=390,y=330)
                baseline.insert(INSERT,columns[1])
                baseline_button=Button(top,image=foto25,bg='#fcfcfc',activebackground='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2')
                baseline_button.place(x=450,y=325)
                baseline_button.bind("<Button-1>",baseline_samples)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',close_window)
        def data_imputation(self):
            del(samples[0:len(samples)])
            def imputation_cancel():
                top.destroy()
                window_status.append('close')
                imputate_list.append('close')
                win.update()
            def replace_method(event):
                def select_method(event):
                    replace.config(state=NORMAL,disabledbackground='#ffffff')
                    replace.delete(0,END)
                    index=replace_list.get(replace_list.curselection())
                    replace.insert(INSERT,index)
                    replace.config(state=DISABLED,disabledbackground='#ffffff')
                    replace_list.destroy()
                    imputate_list.append('close')
                if imputate_list[len(imputate_list)-1]=='close':
                    replace_list=Listbox(top,font=font_type,width=20,height=3, bg='#ffffff',fg='#000000',relief='groove',cursor='hand2')
                    replace_list.place(x=150,y=270)
                    replace_list.insert(END,'Constant','Row average','Column average')
                    replace_list.bind("<Double-Button-1>",select_method)
                    imputate_list.append('open')
            def select_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        selected_columns.insert(END,index)
                        samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_columns(event):
                try:
                    indexes=selected_columns.curselection()
                    for i in range(0,len(indexes),1):
                        index=selected_columns.get(indexes[i])
                        if len(samples)==1:
                            del(samples[0])
                        for i in range (0,len(samples)-1,1):
                            if samples[i]==index:
                                del(samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        selected_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def accept_imputation():
                if len(samples)>=2:
                    if samples.count('GeneID')==0 and samples.count('ID')==0 and samples.count('gene_id')==0:
                        analysis_type.append('imputate')
                        top.destroy()
                        window_status.append('close')
                        top.update()
                        data_table()
                    if samples.count('GeneID')>0 or samples.count('ID')>0 or samples.count('gene_id')>0:
                        warning=showinfo(title='Imputation',message='All columns selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if len(samples)<=1:
                    warning=showinfo(title='Imputation',message='Two or more samples must be selected before data imputation.')
            if window_status[len(window_status)-1]=='close':
                def window_close():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Imputation')
                top.resizable(0,0)
                top.maxsize(520,460)
                Frame(top,bg='#f5f5f5',width=520,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                file_entry.insert(INSERT,file_name[len(file_name)-1])
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                Label(top,text='Selected Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=300,y=40)
                available_columns=Listbox(top,width=30,font=font_type,height=10,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=190,y=60,height=175)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                selected_columns=Listbox(top,width=30,font=font_type,height=10,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                selected_columns.place(x=300,y=60)
                selected_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                selected_columns_scrollbar.place(x=480,y=60,height=175)
                selected_columns.config(yscrollcommand=selected_columns_scrollbar.set)
                selected_columns_scrollbar.config(command=selected_columns.yview)
                available_columns.bind("<Double-Button-1>",select_columns)
                selected_columns.bind("<Double-Button-1>",deselect_columns)
                Label(top,text='Replace empty values for:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=270)
                replace=Entry(top,font=font_type,width=20,bg='#ffffff',fg='#000000',textvariable=replaces,relief='groove')
                replace.place(x=150,y=270)
                replace.config(state=NORMAL,disabledbackground='#ffffff')
                replace.delete(0,END)
                replace.insert(INSERT,'Constant')
                replace.config(state=DISABLED,disabledbackground='#ffffff')
                replace_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                replace_button.place(x=275,y=265)
                replace_button.bind("<Button-1>",replace_method)
                replace_value=Entry(top,font=font_type,width=4,bg='#ffffff',fg='#000000',textvariable=replace_values,relief='groove')
                replace_value.place(x=320,y=270)
                replace_value.delete(0,END)
                replace_value.insert(INSERT,'0')
                accept=Button(top,image=foto26,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                accept.place(x=230,y=110)
                accept.bind("<Button-1>",select_columns)
                accept.bind("<Double-Button-1>chrome",select_columns)
                remove=Button(top,image=foto27,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                remove.place(x=230,y=150)
                remove.bind("<Button-1>",deselect_columns)
                remove.bind("<Double-Button-1>",deselect_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_imputation).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=imputation_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',window_close)
        def data_summarization(self):
            del(samples[0:len(samples)])
            def summarization_cancel():
                top.destroy()
                window_status.append('close')
                win.update()
            def accept_summarization():
                if len(samples)>=2:
                    if samples.count('GeneID')==0 and samples.count('ID')==0 and samples.count('gene_id')==0:
                        analysis_type.append('summarize')
                        top.destroy()
                        window_status.append('close')
                        top.update()
                        data_table()
                    if samples.count('GeneID')>0 or samples.count('ID')>0 or samples.count('gene_id')>0:
                        warning=showinfo(title='Summarization',message='All columns selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if len(samples)<=1:
                    warning=showinfo(title='Summarization',message='Two or more samples must be selected before data summarization.')
            def replace_method(event):
                def select_method(event):
                    replace.delete(0,END)
                    index=replace_list.get(replace_list.curselection())
                    replace.insert(INSERT,index)
                    replace_list.destroy()
                replace_list=Listbox(top,font=font_type,width=20,height=4,bg='#fcfcfc',fg='#000000',relief='groove',cursor='hand2')
                replace_list.place(x=150,y=270)
                replace_list.insert(END,'Average','Median','Standard deviation','Variance')
                replace_list.bind("<Double-Button-1>",select_method)
            def select_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        selected_columns.insert(END,index)
                        samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_columns(event):
                try:
                    indexes=selected_columns.curselection()
                    for i in range(0,len(indexes),1):
                        index=selected_columns.get(indexes[i])
                        if len(samples)==1:
                            del(samples[0])
                        for i in range (0,len(samples)-1,1):
                            if samples[i]==index:
                                del(samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        selected_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Summarize Data')
                top.resizable(0,0)
                top.maxsize(520,460)
                Frame(top,bg='#f5f5f5',width=520,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                file_entry.insert(INSERT,file_name[len(file_name)-1])
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                Label(top,text='Selected Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=300,y=40)
                available_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=190,y=60,height=175)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                selected_columns=Listbox(top,width=30,font=font_type,height=10,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                selected_columns.place(x=300,y=60)
                selected_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                selected_columns_scrollbar.place(x=480,y=60,height=175)
                selected_columns.config(yscrollcommand=selected_columns_scrollbar.set)
                selected_columns_scrollbar.config(command=selected_columns.yview)
                available_columns.bind("<Double-Button-1>",select_columns)
                selected_columns.bind("<Double-Button-1>",deselect_columns)
                Label(top,text='Summarization measure:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=270)
                replace=Entry(top,font=font_type,width=20,bg='#fcfcfc',fg='#000000',textvariable=summarization_replaces,relief='groove')
                replace.place(x=150,y=270)
                replace.delete(0,END)
                replace.insert(INSERT,'Average')
                replace_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                replace_button.place(x=275,y=265)
                replace_button.bind("<Button-1>",replace_method)
                accept=Button(top,image=foto26,bg='#f5f5f5',relief='flat',activebackground='#f5f5f5',overrelief='groove',cursor='hand2')
                accept.place(x=230,y=110)
                accept.bind("<Button-1>",select_columns)
                accept.bind("<Double-Button-1>",select_columns)
                remove=Button(top,image=foto27,bg='#f5f5f5',activebackground='#f5f5f5',cursor='hand2',relief='flat',overrelief='groove')
                remove.place(x=230,y=150)
                remove.bind("<Button-1>",deselect_columns)
                remove.bind("<Double-Button-1>",deselect_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_summarization).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=summarization_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',close_window)
        def data_transposition(self):
            if genomic_technology[len(genomic_technology)-1]=='Custom':
                if len(data1)<=100:
                    del(gene_id[0:len(gene_id)]);del(data1[0:len(data1)]);del(data2[0:len(data2)])
                    del(new_data1[0:len(new_data1)]);del(new_data2[0:len(new_data2)]);del(values[0:len(values)]);del(values2[0:len(values2)])
                    values.append(0),values2.append(0)
                    infiles=open('c:\\PG\\processes\\data\\samples.txt','w')
                    for i in range(1,len(columns),1):
                        infiles.write(columns[i]+'\n')
                    infiles.close()
                    r.source('C:\\PG\\srs.R')
                    r.transpose('C:\\PG\\processes\\data\\result.txt')
                    file='C:\\PG\\processes\\data\\result.txt'
                    content=open(file,'r')
                    contents=content.readlines()
                    column_number=(contents[0].count('\t'))+1
                    infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
                    infile=r.matrix(infile)
                    del (columns[0:len(columns)])
                    for i in range(0,len(infile),1):
                        columns.append(infile[i][0][0])
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                    win.update()
                    for i in range(0,len(data1),1):
                        if data1[i]>values[len(values)-1]:
                            values.append(data1[i])
                        if data2[i]>values2[len(values2)-1]:
                            values2.append(data2[i])
                        if values2[len(values2)-1]>values[len(values)-1]:
                            values.append(values2[len(values2)-1])
                    for i in range(0,len(data1),1):
                        new_data1.append((data1[i]*650)/values[len(values)-1])
                        new_data2.append((data2[i]*400)/values[len(values)-1])
                    analysis_type.append('raw')
                    data_transpose.append('yes')
                    data_table()
                    win.update()
                    sample1.delete(0,END)
                    sample2.delete(0,END)
                    sample1.insert(INSERT,columns[2])
                    sample2.insert(INSERT,columns[1])
                    graphic()
            if len(data1)>100:
                warning=showinfo(title='Data Transposition',message='Table size makes impossible to transpose')
                data_transpose.append('no')
    # METODOS ESTADISTICOS DE EXPRESION DIFERENCIAL
    class statistic_class():
        def segmentation_process(self):
            if analysis_type.count('segmentation')==0 and analysis_type.count('tidy_up')>=1:
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.resizable(0,0)
                top.title('Segmentation')
                def segmentation_methods():
                      def select_segmentation(event):
                            index=method_list.get(method_list.curselection())
                            method_list.destroy()
                            segmentation_method.delete(0,END)
                            segmentation_method.insert(INSERT,index)
                      method_list=Listbox(top,width=20,height=2,font=font_type,bg='#ffffff',relief='groove',cursor='hand2')
                      method_list.place(x=200,y=350)
                      method_list.insert(END,'DNAcopy','GLAD')
                      method_list.bind("<Double-Button-1>",select_segmentation)
                def close_window():
                    top.destroy()
                    window_status.append('close')
                def select_control_columns(event):
                    try:
                        indexes=available_columns.curselection()
                        for i in range (0,len(indexes),1):
                            index=available_columns.get(indexes[i])
                            control_columns.insert(END,index)
                            control_samples.append(index)
                        for i in range (0,len(indexes),1):
                            available_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                def deselect_control_columns(event):
                    try:
                        indexes=control_columns.curselection()
                        for j in range(0,len(indexes),1):
                            index=control_columns.get(indexes[j])
                            if len(control_samples)==1:
                                del(control_samples[0])
                            for i in range (0,len(control_samples)-1,1):
                                if control_samples[i]==index:
                                    del(control_samples[i])
                            available_columns.insert(END,index)
                        for i in range(0,len(indexes),1):
                            control_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                def select_target_columns(event):
                    try:
                        indexes=available_columns.curselection()
                        for i in range (0,len(indexes),1):
                            index=available_columns.get(indexes[i])
                            target_columns.insert(END,index)
                            target_samples.append(index)
                        for i in range (0,len(indexes),1):
                            available_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                def deselect_target_columns(event):
                    try:
                        indexes=target_columns.curselection()
                        for j in range(0,len(indexes),1):
                            index=target_columns.get(indexes[j])
                            if len(target_samples)==1:
                                del(target_samples[0])
                            for i in range (0,len(target_samples)-1,1):
                                if target_samples[i]==index:
                                    del(target_samples[i])
                            available_columns.insert(END,index)
                        for i in range(0,len(indexes),1):
                            target_columns.delete(int(indexes[i])-i)
                    except:
                        warning=showwarning(title='Select Columns',message='No column selected')
                def accept_segmentation():
                    if control_samples.count('GeneID')==0 and control_samples.count('ID')==0 and control_samples.count('gene_id')==0 and target_samples.count('GeneID')==0 and target_samples.count('ID')==0 and target_samples.count('gene_id')==0:
                        segmentation_list.append(segmentation_method.get())
                        workflow.append('segmentation')
                        analysis_type.append('segmentation')
                        for i in range(0,len(control_samples),1):
                            samples.append(control_samples[i])
                            design.append(1)
                        for i in range(0,len(target_samples),1):
                            samples.append(target_samples[i])
                            design.append(-1)
                        infile=open('c:\\PG\\processes\\data\\design.txt','w')
                        for i in range(0,len(design),1):
                            infile.write(str(design[i])+'\n')
                        infile.close()
                        top.destroy()
                        window_status.append('close')
                        win.update()
                        data_table()
                    if control_samples.count('GeneID')>0 or control_samples.count('ID')>0 or control_samples.count('gene_id')>0 or target_samples.count('GeneID')>0 or target_samples.count('ID')>0 or target_samples.count('gene_id')>0:
                        warning=showinfo(title='Segmentation',message='All samples selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                Frame(top,width=480,height=500,bg='#f5f5f5',borderwidth=2).pack()
                Label(top,text='Segmentation process divides the genome according the estimation of significant alterations',font=font_type,bg='#f5f5f5').place(x=10,y=10)
                Label(top,text='Samples',font=font_type,bg='#f5f5f5').place(x=10,y=60)
                Label(top,text='Control',font=font_type,bg='#f5f5f5').place(x=10,y=200)
                Label(top,text='Target',font=font_type,bg='#f5f5f5').place(x=300,y=200)
                Label(top,text='Select the Segmentation method:',font=font_type,bg='#f5f5f5').place(x=10,y=350)
                Label(top,text='Select the Gain and Loss threshold:',font=font_type,bg='#f5f5f5').place(x=10,y=390)
                available_columns=Listbox(top,width=71,height=6,font=font_type,relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=80)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=440,y=80,height=105)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                control_columns=Listbox(top,width=25,height=6,font=font_type,relief='groove',selectmode=MULTIPLE)
                control_columns.place(x=10,y=220)
                target_columns=Listbox(top,width=25,height=6,font=font_type,relief='groove',selectmode=MULTIPLE)
                target_columns.place(x=300,y=220)
                segmentation_method=Entry(top,width=20,font=font_type,relief='groove')
                segmentation_method.place(x=200,y=350)
                segmentation_method.delete(0,END)
                segmentation_method.insert(INSERT,'DNAcopy')
                gain_threshold=Entry(top,width=20,font=font_type,relief='groove')
                gain_threshold.place(x=200,y=390)
                gain_threshold.config(state=DISABLED,disabledbackground='#f5f5f5')
                deselect_control=Button(top,image=foto24,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_control.place(x=100,y=190)
                deselect_control.bind("<Button-1>",deselect_control_columns)
                select_control=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_control.place(x=130,y=190)
                select_control.bind("<Button-1>",select_control_columns)
                deselect_target=Button(top,image=foto24,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_target.place(x=390,y=190)
                deselect_target.bind("<Button-1>",deselect_target_columns)
                select_target=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_target.place(x=420,y=190)
                select_target.bind("<Button-1>",select_target_columns)
                Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=segmentation_methods).place(x=325,y=345)
                Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2').place(x=325,y=385)
                Button(top,text='Help',font=font_type,width=7,bg='#ffffff',relief='groove',cursor='hand2').place(x=10,y=460)
                Button(top,text='OK',font=font_type,width=7,bg='#ffffff',relief='groove',cursor='hand2',command=accept_segmentation).place(x=360,y=460)
                Button(top,text='Cancel',font=font_type,width=7,bg='#ffffff',relief='groove',cursor='hand2',command=close_window).place(x=420,y=460)
                for i in range(0,len(columns),1):
                        available_columns.insert(END,columns[i])
                top.protocol('WM_DELETE_WINDOW',close_window)               
            if analysis_type.count('segmentation')>=1:
                warning=showwarning(title='Segmentation',message='Segmentation process has already been applied')
        def limma(self):
            del(samples[0:len(samples)])
            del(control_samples[0:len(control_samples)])
            del(target_samples[0:len(target_samples)])
            del(design[0:len(design)])
            def limma_cancel():
                top.destroy()
                window_status.append('close')
                win.update()
            def select_control_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        control_columns.insert(END,index)
                        control_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_control_columns(event):
                try:
                    indexes=control_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=control_columns.get(indexes[j])
                        if len(control_samples)==1:
                            del(control_samples[0])
                        for i in range (0,len(control_samples)-1,1):
                            if control_samples[i]==index:
                                del(control_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        control_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def select_target_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        target_columns.insert(END,index)
                        target_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_target_columns(event):
                try:
                    indexes=target_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=target_columns.get(indexes[j])
                        if len(target_samples)==1:
                            del(target_samples[0])
                        for i in range (0,len(target_samples)-1,1):
                            if target_samples[i]==index:
                                del(target_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        target_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def accept_limma():
                if len(control_samples)>=2 and len(target_samples)>=2 or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    if control_samples.count('GeneID')==0 and control_samples.count('ID')==0 and control_samples.count('gene_id')==0 and target_samples.count('GeneID')==0 and target_samples.count('ID')==0 and target_samples.count('gene_id')==0:
                        analysis_type.append('limma')
                        if len(control_samples)>=2:
                            for i in range(0,len(control_samples),1):
                                samples.append(control_samples[i])
                                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                                    design.append(0)
                                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                                    design.append(-1)
                        if len(target_samples)>=2:
                            for i in range(0,len(target_samples),1):
                                samples.append(target_samples[i])
                                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                                    design.append(1)
                                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                                    design.append(0)
                        if len(control_samples)>=2 or len(target_samples)>=2:
                            infile=open('c:\\PG\\processes\\data\\design.txt','w')
                            for i in range(0,len(design),1):
                                infile.write(str(design[i])+'\n')
                            infile.close()
                            top.destroy()
                            window_status.append('close')
                            win.update()
                            data_table()
                    if control_samples.count('GeneID')>0 or control_samples.count('ID')>0 or control_samples.count('gene_id')>0 or target_samples.count('GeneID')>0 or target_samples.count('ID')>0 or target_samples.count('gene_id')>0:
                        warning=showinfo(title='LIMMA',message='All samples selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' and genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                    if len(control_samples)<=1 or len(target_samples)<=1:
                        warning=showinfo(title='LIMMA',message='Two or more samples must be selected as target and control before applying LIMMA analysis.')
                        window_status.append('close')
                        top.destroy()
            def fdr_methods(event):
                def fdr_selection(event):
                    method_selected=fdr_list.get(fdr_list.curselection())
                    fdr_list.destroy()
                    fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                    fdr_method.delete(0,END)
                    fdr_method.insert(INSERT,method_selected)
                    fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                    fdr_method_list.append('close')
                if fdr_method_list[len(fdr_method_list)-1]=='close':
                    fdr_method_list.append('open')
                    fdr_list=Listbox(top,width=30,height=4,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove')
                    fdr_list.place(x=180,y=330)
                    fdr_list.insert(END,'Benjamini-Hochberg','Holm','Bonferroni','none')
                    fdr_list.bind("<Double-Button-1>",fdr_selection)
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Linear Models for Microarrays')
                top.resizable(0,0)
                top.maxsize(500,460)
                Frame(top,bg='#f5f5f5',width=500,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                file_entry.config(state=NORMAL,disabledbackground='#f5f5f5')
                file_entry.insert(INSERT,file_name[len(file_name)-1])
                file_entry.config(state=DISABLED,disabledbackground='#f5f5f5')
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                available_columns=Listbox(top,width=77,font=font_type,height=7,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=475,y=60,height=120)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                Label(top,text='Control',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=200)
                Label(top,text='Target',font=font_type,bg='#f5f5f5',relief='flat').place(x=260,y=200)
                control_columns=Listbox(top,width=35,font=font_type,height=4,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                control_columns.place(x=10,y=220)
                target_columns=Listbox(top,width=35,font=font_type,height=4,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                target_columns.place(x=260,y=220)
                select_control=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_control.place(x=100,y=190)
                select_control.bind("<Button-1>",select_control_columns)
                deselect_control=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_control.place(x=130,y=190)
                deselect_control.bind("<Button-1>",deselect_control_columns)
                select_target=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_target.place(x=350,y=190)
                select_target.bind("<Button-1>",select_target_columns)
                deselect_target=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_target.place(x=380,y=190)
                deselect_target.bind("<Button-1>",deselect_target_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_limma).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=limma_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                Label(top,text='False Discovery Rate method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=330)
                fdr_method=Entry(top,font=font_type,bg='#ffffff',width=30,textvariable=fdr,relief='groove')
                fdr_method.place(x=180,y=330)
                fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                fdr_method.delete(0,END)
                fdr_method.insert(INSERT,'Benjamini-Hochberg')
                fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                fdr_method_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                fdr_method_button.place(x=365,y=325)
                fdr_method_button.bind("<Button-1>",fdr_methods)
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',close_window)
        def t_test(self):
            del(samples[0:len(samples)])
            del(control_samples[0:len(control_samples)])
            del(target_samples[0:len(target_samples)])
            def ttest_cancel():
                top.destroy()
                window_status.append('close')
                win.update()
            def load_design_file():
                file=askopenfilename()
                design_file.delete(0,END)
                design_file(INSERT,file)
            def select_control_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        control_columns.insert(END,index)
                        control_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_control_columns(event):
                try:
                    indexes=control_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=control_columns.get(indexes[j])
                        if len(control_samples)==1:
                            del(control_samples[0])
                        for i in range (0,len(control_samples)-1,1):
                            if control_samples[i]==index:
                                del(control_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        control_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def select_target_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        target_columns.insert(END,index)
                        target_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_target_columns(event):
                try:
                    indexes=target_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=target_columns.get(indexes[j])
                        if len(target_samples)==1:
                            del(target_samples[0])
                        for i in range (0,len(target_samples)-1,1):
                            if target_samples[i]==index:
                                del(target_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        target_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def fdr_methods(event):
                def fdr_selection(event):
                    method_selected=fdr_list.get(fdr_list.curselection())
                    fdr_list.destroy()
                    fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                    fdr_method.delete(0,END)
                    fdr_method.insert(INSERT,method_selected)
                    fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                    fdr_method_list.append('close')
                if fdr_method_list[len(fdr_method_list)-1]=='close':
                    fdr_method_list.append('open')
                    fdr_list=Listbox(top,width=30,height=4,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove')
                    fdr_list.place(x=180,y=360)
                    fdr_list.insert(END,'Benjamini-Hochberg','Holm','Bonferroni','none')
                    fdr_list.bind("<Double-Button-1>",fdr_selection)
            def accept_ttest():
                if len(control_samples)>=2 and len(target_samples)>=2 and genomic_technology[len(genomic_technology)-1]!='Genepix' and genomic_technology[len(genomic_technology)-1]!='Quantarray' and genomic_technology[len(genomic_technology)-1]!='Scanarray' and genomic_technology[len(genomic_technology)-1]!='Agilent' and genomic_technology[len(genomic_technology)-1]!='Imagene':
                    if control_samples.count('GeneID')==0 and control_samples.count('ID')==0 and control_samples.count('gene_id')==0 and target_samples.count('GeneID')==0 and target_samples.count('ID')==0 and target_samples.count('gene_id')==0:
                        analysis_type.append('ttest')
                        for i in range(0,len(control_samples),1):
                            samples.append(control_samples[i])
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina':
                                design.append(0)
                        for i in range(0,len(target_samples),1):
                            samples.append(target_samples[i])
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina':
                                design.append(1)
                        infile=open('c:\\PG\\processes\\data\\design.txt','w')
                        for i in range(0,len(design),1):
                            infile.write(str(design[i])+'\n')
                        infile.close()
                        top.destroy()
                        window_status.append('close')
                        win.update()
                        data_table()
                    if control_samples.count('GeneID')>0 or control_samples.count('ID')>0 or control_samples.count('gene_id')>0 or target_samples.count('GeneID')>0 or target_samples.count('ID')>0 or target_samples.count('gene_id')>0:
                        warning=showinfo(title='T-Test',message='All columns selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if len(control_samples)<=1 or len(target_samples)<=1:
                    warning=showinfo(title='T-Test',message='Two or more samples must be selected as target and control before applying T-Test analysis.')
                    window_status.append('close')
                    top.destroy()
                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    warning=showinfo(title='T-Test',message='Two channel technologies do not support T-Test analysis.')
                    window_status.append('close')
                    top.destroy()
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Unpaired T-Test Analysis')
                top.resizable(0,0)
                top.maxsize(500,460)
                Frame(top,bg='#f5f5f5',width=500,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                file_entry.config(state=NORMAL,disabledbackground='#f5f5f5')
                file_entry.insert(INSERT,file_name[len(file_name)-1])
                file_entry.config(state=DISABLED,disabledbackground='#f5f5f5')
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                available_columns=Listbox(top,width=77,font=font_type,height=7,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=475,y=60,height=120)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                Label(top,text='Control',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=200)
                Label(top,text='Target',font=font_type,bg='#f5f5f5',relief='flat').place(x=260,y=200)
                control_columns=Listbox(top,width=35,font=font_type,height=4,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                control_columns.place(x=10,y=220)
                target_columns=Listbox(top,width=35,font=font_type,height=4,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                target_columns.place(x=260,y=220)
                select_control=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_control.place(x=100,y=190)
                select_control.bind("<Button-1>",select_control_columns)
                deselect_control=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_control.place(x=130,y=190)
                deselect_control.bind("<Button-1>",deselect_control_columns)
                select_target=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_target.place(x=350,y=190)
                select_target.bind("<Button-1>",select_target_columns)
                deselect_target=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_target.place(x=380,y=190)
                deselect_target.bind("<Button-1>",deselect_target_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_ttest).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=ttest_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                Label(top,text='False Discovery Rate method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=360)
                fdr_method=Entry(top,font=font_type,bg='#ffffff',width=30,textvariable=fdr,relief='groove')
                fdr_method.place(x=180,y=360)
                fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                fdr_method.delete(0,END)
                fdr_method.insert(INSERT,'Benjamini-Hochberg')
                fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                fdr_method_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                fdr_method_button.place(x=365,y=355)
                fdr_method_button.bind("<Button-1>",fdr_methods)
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',close_window)
        def sam(self):
            del(samples[0:len(samples)])
            del(control_samples[0:len(control_samples)])
            del(target_samples[0:len(target_samples)])
            def sam_cancel():
                top.destroy()
                window_status.append('close')
                win.update()
            def load_design_file():
                file=askopenfilename()
                design_file.delete(0,END)
                design_file(INSERT,file)
            def select_control_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        control_columns.insert(END,index)
                        control_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_control_columns(event):
                try:
                    indexes=control_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=control_columns.get(indexes[j])
                        if len(control_samples)==1:
                            del(control_samples[0])
                        for i in range (0,len(control_samples)-1,1):
                            if control_samples[i]==index:
                                del(control_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        control_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def select_target_columns(event):
                try:
                    indexes=available_columns.curselection()
                    for i in range (0,len(indexes),1):
                        index=available_columns.get(indexes[i])
                        target_columns.insert(END,index)
                        target_samples.append(index)
                    for i in range (0,len(indexes),1):
                        available_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def deselect_target_columns(event):
                try:
                    indexes=target_columns.curselection()
                    for j in range(0,len(indexes),1):
                        index=target_columns.get(indexes[j])
                        if len(target_samples)==1:
                            del(target_samples[0])
                        for i in range (0,len(target_samples)-1,1):
                            if target_samples[i]==index:
                                del(target_samples[i])
                        available_columns.insert(END,index)
                    for i in range(0,len(indexes),1):
                        target_columns.delete(int(indexes[i])-i)
                except:
                    warning=showwarning(title='Select Columns',message='No column selected')
            def fdr_methods(event):
                def fdr_selection(event):
                    method_selected=fdr_list.get(fdr_list.curselection())
                    fdr_list.destroy()
                    fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                    fdr_method.delete(0,END)
                    fdr_method.insert(INSERT,method_selected)
                    fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                    fdr_method_list.append('close')
                if fdr_method_list[len(fdr_method_list)-1]=='close':
                    fdr_method_list.append('open')
                    fdr_list=Listbox(top,width=30,height=4,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove')
                    fdr_list.place(x=180,y=360)
                    fdr_list.insert(END,'Benjamini-Hochberg','Holm','Bonferroni','none')
                    fdr_list.bind("<Double-Button-1>",fdr_selection)
            def statistic_methods(event):
                def statistic_selection(event):
                    method_selected=statistic_list.get(statistic_list.curselection())
                    statistic_list.destroy()
                    statistic_method.config(state=NORMAL,disabledbackground='#ffffff')
                    statistic_method.delete(0,END)
                    statistic_method.insert(INSERT,method_selected)
                    statistic_method.config(state=DISABLED,disabledbackground='#ffffff')
                    sam_statistic.append('close')
                if sam_statistic[len(sam_statistic)-1]=='close':
                    sam_statistic.append('open')
                    statistic_list=Listbox(top,width=30,height=2,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove')
                    statistic_list.place(x=180,y=330)
                    statistic_list.insert(END,'Standard','Wilcoxon')
                    statistic_list.bind("<Double-Button-1>",statistic_selection)
            def accept_sam():
                if len(control_samples)>=2 and len(target_samples)>=2 and genomic_technology[len(genomic_technology)-1]!='Genepix' and genomic_technology[len(genomic_technology)-1]!='Quantarray' and genomic_technology[len(genomic_technology)-1]!='Scanarray' and genomic_technology[len(genomic_technology)-1]!='Agilent' and genomic_technology[len(genomic_technology)-1]!='Imagene':
                    if control_samples.count('GeneID')==0 and control_samples.count('ID')==0 and control_samples.count('gene_id')==0 and target_samples.count('GeneID')==0 and target_samples.count('ID')==0 and target_samples.count('gene_id')==0:
                        analysis_type.append('sam')
                        for i in range(0,len(control_samples),1):
                            samples.append(control_samples[i])
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina':
                                design.append(0)
                        for i in range(0,len(target_samples),1):
                            samples.append(target_samples[i])
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina':
                                design.append(1)
                        infile=open('c:\\PG\\processes\\data\\design.txt','w')
                        for i in range(0,len(design),1):
                            infile.write(str(design[i])+'\n')
                        infile.close()
                        top.destroy()
                        window_status.append('close')
                        win.update()
                        data_table()
                    if control_samples.count('GeneID')>0 or control_samples.count('ID')>0 or control_samples.count('gene_id')>0 or target_samples.count('GeneID')>0 or target_samples.count('ID')>0 or target_samples.count('gene_id')>0:
                        warning=showinfo(title='SAM',message='All columns selected must be numeric.')
                        window_status.append('close')
                        top.destroy()
                if len(control_samples)<=1 or len(target_samples)<=1:
                    warning=showinfo(title='SAM',message='Two or more samples must be selected as target and control before applying SAM analysis.')
                    window_status.append('close')
                    top.destroy()
                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    warning=showinfo(title='SAM',message='Two channel technologies do not support SAM analysis.')
                    window_status.append('close')
                    top.destroy()
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.title('Significant Analysis of Microarrays')
                top.resizable(0,0)
                top.maxsize(500,460)
                Frame(top,bg='#f5f5f5',width=500,height=460,relief='groove').pack()
                Label(top,text='File Selected:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=10)
                file_entry=Entry(top,font=font_type,width=40,bg='#f5f5f5',fg='#0000ff',relief='flat')
                file_entry.place(x=80,y=10)
                file_entry.config(state=NORMAL,disabledbackground='#f5f5f5')
                file_entry.insert(INSERT,file_name[len(file_name)-1])
                file_entry.config(state=DISABLED,disabledbackground='#f5f5f5')
                Label(top,text='Data Columns',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=40)
                available_columns=Listbox(top,width=77,font=font_type,height=7,bg='#ffffff',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                available_columns.place(x=10,y=60)
                available_columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                available_columns_scrollbar.place(x=475,y=60,height=120)
                available_columns.config(yscrollcommand=available_columns_scrollbar.set)
                available_columns_scrollbar.config(command=available_columns.yview)
                Label(top,text='Control',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=200)
                Label(top,text='Target',font=font_type,bg='#f5f5f5',relief='flat').place(x=260,y=200)
                control_columns=Listbox(top,width=35,font=font_type,height=4,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                control_columns.place(x=10,y=220)
                target_columns=Listbox(top,width=35,font=font_type,height=4,bg='#fcfcfc',cursor='hand2',relief='groove',selectmode=MULTIPLE)
                target_columns.place(x=260,y=220)
                select_control=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_control.place(x=100,y=190)
                select_control.bind("<Button-1>",select_control_columns)
                deselect_control=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_control.place(x=130,y=190)
                deselect_control.bind("<Button-1>",deselect_control_columns)
                select_target=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                select_target.place(x=350,y=190)
                select_target.bind("<Button-1>",select_target_columns)
                deselect_target=Button(top,image=foto24,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                deselect_target.place(x=380,y=190)
                deselect_target.bind("<Button-1>",deselect_target_columns)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_sam).place(x=370,y=420)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=sam_cancel).place(x=425,y=420)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=420)
                for i in range(0,len(columns),1):
                    available_columns.insert(END,columns[i])
                Label(top,text='False Discovery Rate method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=360)
                fdr_method=Entry(top,font=font_type,bg='#ffffff',width=30,textvariable=fdr,relief='groove')
                fdr_method.place(x=180,y=360)
                fdr_method.config(state=NORMAL,disabledbackground='#ffffff')
                fdr_method.delete(0,END)
                fdr_method.insert(INSERT,'Benjamini-Hochberg')
                fdr_method.config(state=DISABLED,disabledbackground='#ffffff')
                fdr_method_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                fdr_method_button.place(x=365,y=355)
                fdr_method_button.bind("<Button-1>",fdr_methods)
                Label(top,text='Statistic method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=10,y=330)
                statistic_method=Entry(top,font=font_type,bg='#ffffff',width=30,textvariable=statistic,relief='groove')
                statistic_method.place(x=180,y=330)
                statistic_method.config(state=NORMAL,disabledbackground='#ffffff')
                statistic_method.delete(0,END)
                statistic_method.insert(INSERT,'Standard')
                statistic_method.config(state=DISABLED,disabledbackground='#ffffff')
                statistic_method_button=Button(top,image=foto25,bg='#f5f5f5',activebackground='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                statistic_method_button.place(x=365,y=325)
                statistic_method_button.bind("<Button-1>",statistic_methods)
                window_status.append('open')
                top.protocol('WM_DELETE_WINDOW',close_window)
    #MODELOS DE CARGA DE TECNOLOGIAS SOPORTADAS
    class technology_class():
        def affymetrix(self):
            def summarize_affymetrix(event):
                if workflow.count('preprocessing')>0:
                    warning=showinfo(title='Workflow',message='Microarrays already preprocessed.')
                if workflow.count('preprocessing')==0:
                    workflow.append('preprocessing')
                    if analysis_method.get()==1:
                        affymetrix_method.append('RMA')
                        analysis_type.append('All-Affymetrix')
                        genomic_technology.append('All-Affymetrix')
                    if analysis_method.get()==0:
                        affymetrix_method.append('MAS5')
                    data_table()
                    information_text.config(state=NORMAL)
                    if analysis_method.get()==1:
                        information_text.insert(INSERT,'- RMA method used to preprocess the Affymetrix data. Check quality controls to decide if normalization is required.\n')
                    if analysis_method.get()==0:
                        information_text.insert(INSERT,'- MAS5.0 method used to preprocess the Affymetrix data. Normalization step is mandatory.\n')
                    information_text.config(state=DISABLED)
            def affymetrix_normalization(event):
                if workflow.count('preprocessing')>=1:
                    a=preprocess_class()
                    a.data_normalization()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before normalizing')
                workflow.append('normalization')
            def affymetrix_limma(event):
                if workflow.count('preprocessing')>=1:
                    workflow.append('limma')
                    a=statistic_class()
                    a.limma()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before differential expression analysis')
            def affymetrix_quality(event):
                if workflow.count('preprocessing')>=1:
                    a=quality_class()
                    a.array_quality()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before quality visualizations')
            def affymetrix_annotation(event):
                if workflow.count('limma')>=1:
                    a=functional_analysis()
                    a.retrieve_annotations()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Differential expression must be applied before gene annotation')
            def affymetrix_clustering(event):
                if workflow.count('limma')>=1:
                    a=preprocess_class()
                    a.hierarchical_clustering()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
            presentation_text.destroy()  
            Label(win,text='AFFYMETRIX DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)            
            preprocess=Label(win,text='1. Microarrays Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            preprocess.place(x=7,y=100)
            preprocess.bind("<Button-1>",summarize_affymetrix)
            quality=Label(win,text='2. Quality Control',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            quality.place(x=7,y=130)
            quality.bind("<Button-1>",affymetrix_quality)
            normalization=Label(win,text='3. Normalization',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            normalization.place(x=7,y=160)
            normalization.bind("<Button-1>",affymetrix_normalization)
            differential_expression=Label(win,text='4. Expression Analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            differential_expression.place(x=7,y=190)
            differential_expression.bind("<Button-1>",affymetrix_limma)
            functional=Label(win,text='5. Annotation Retriever',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            functional.place(x=7,y=220)
            functional.bind("<Button-1>",affymetrix_annotation)
            cluster=Label(win,text='6. Hierarchical Clustering',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            cluster.place(x=7,y=250)
            cluster.bind("<Button-1>",affymetrix_clustering)
            win.update()
        def custom (self):
            def custom_normalization(event):
                a=preprocess_class()
                a.data_normalization()
            def custom_limma(event):
                a=statistic_class()
                a.limma()
            def custom_quality(event):
                a=quality_class()
                a.array_quality()
            def custom_annotation(event):
                if workflow.count('limma')>=1:
                    a=functional_analysis()
                    a.retrieve_annotations()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Differential expression must be applied before gene annotation')
            def custom_clustering(event):
                if workflow.count('limma')>=1:
                    a=preprocess_class()
                    a.hierarchical_clustering()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
            presentation_text.destroy()
            Label(win,text='CUSTOM DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            quality=Label(win,text='1. Quality Control',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            quality.place(x=7,y=100)
            quality.bind("<Button-1>",custom_quality)
            normalization=Label(win,text='2. Normalization',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            normalization.place(x=7,y=130)
            normalization.bind("<Button-1>",custom_normalization)
            differential_expression=Label(win,text='3. Expression Analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            differential_expression.place(x=7,y=160)
            differential_expression.bind("<Button-1>", custom_limma)
            functional=Label(win,text='4. Annotation Retriever',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            functional.place(x=7,y=190)
            functional.bind("<Button-1>",custom_annotation)
            Cluster=Label(win,text='5. Hierarchical Clustering',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            Cluster.place(x=7,y=220)
            Cluster.bind("<Button-1>", custom_clustering)
            win.update()
        def illumina(self):
            presentation_text.destroy()  
            def summarize_illumina(event):
                if workflow.count('preprocessing')>0:
                    warning=showinfo(title='Workflow',message='Microarrays already preprocessed.')
                if workflow.count('preprocessing')==0:
                    workflow.append('preprocessing')
                    data_table()
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- Illumina raw data imported. Normalization is highly recommended.\n')
                    information_text.config(state=DISABLED)
            def illumina_quality(event):
                if workflow.count('preprocessing')>=1:
                    a=quality_class()
                    a.array_quality()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before quality visualizations')
            def illumina_normalization(event):
                if workflow.count('preprocessing')>=1:
                    a=preprocess_class()
                    a.data_normalization()
                    workflow.append('normalization')
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before quality visualizations')
            def illumina_limma(event):
                if workflow.count('preprocessing')>=1:
                    workflow.append('limma')
                    a=statistic_class()
                    a.limma()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before quality visualizations')
            def illumina_annotation(event):
                if workflow.count('limma')>=1:
                    a=functional_analysis()
                    a.retrieve_annotations()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Differential expression must be applied before gene annotation')
            def illumina_clustering(event):
                if workflow.count('limma')>=1:
                    a=preprocess_class()
                    a.hierarchical_clustering()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
            Label(win,text='ILLUMINA DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            preprocess=Label(win,text='1. Microarrays Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            preprocess.place(x=7,y=100)
            preprocess.bind("<Button-1>",summarize_illumina)
            quality=Label(win,text='2. Quality Control',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            quality.place(x=7,y=130)
            quality.bind("<Button-1>",illumina_quality)
            normalization=Label(win,text='3. Normalization',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            normalization.place(x=7,y=160)
            normalization.bind("<Button-1>",illumina_normalization)
            differential_expression=Label(win,text='4. Expression Analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            differential_expression.place(x=7,y=190)
            differential_expression.bind("<Button-1>",illumina_limma)
            functional=Label(win,text='5. Annotation Retriever',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            functional.place(x=7,y=220)
            functional.bind("<Button-1>",illumina_annotation)
            Cluster=Label(win,text='6. Hierarchical Clustering',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            Cluster.place(x=7,y=250)
            Cluster.bind("<Button-1>", illumina_clustering)
            win.update()
        def two_channel(self):
            presentation_text.destroy()
            if tidy_up_mark[len(tidy_up_mark)-1]=='TRUE':
                def genome_view(event):
                    r.source('c:\\PG\\srs.R')
                    r.tidyview('c:\\PG\\processes\\data\\result.txt')
                tidy_up_view_2=Label(win,text='view',font=font_type,bg='#fcfcfc',fg='#8080c0',relief='flat',cursor='hand2')
                tidy_up_view_2.place(x=120,y=195)
                tidy_up_view_2.bind("<Button-1>",genome_view)
            def summarize_two_channel(event):
                if workflow.count('preprocessing')>0:
                    warning=showinfo(title='Workflow',message='Microarrays already preprocessed.')
                if workflow.count('preprocessing')==0:
                    workflow.append('preprocessing')
                    data_table()
                    information_text.config(state=NORMAL)
                    if genomic_technology[len(genomic_technology)-1]=='Genepix':
                        information_text.insert(INSERT,'- Genepix data imported. Normalization is highly recommended.\n')
                    if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                        information_text.insert(INSERT,'- Agilent data imported. Normalization is highly recommended.\n')
                    if genomic_technology[len(genomic_technology)-1]=='Imagene':
                        information_text.insert(INSERT,'- Imagene data imported. Normalization is highly recommended.\n')
                    if genomic_technology[len(genomic_technology)-1]=='Scanarray':
                        information_text.insert(INSERT,'- Scanarray data imported. Normalization is highly recommended.\n')
                    if genomic_technology[len(genomic_technology)-1]=='Quantarray':
                        information_text.insert(INSERT,'- Quantarray data imported. Normalization is highly recommended.\n')
                    information_text.config(state=DISABLED)
            def two_channel_quality(event):
                if workflow.count('preprocessing')>=1:
                    a=quality_class()
                    a.array_quality()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Microarrays must be preprocessed before quality visualizations')
            def two_channel_normalization(event):
                if workflow.count('preprocessing')>=1:
                    a=preprocess_class()
                    a.data_normalization()
                    workflow.append('normalization')
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Normalization',message='Microarrays must be preprocessed before quality visualizations')
            def two_channel_limma(event):
                if workflow.count('preprocessing')>=1:
                    workflow.append('limma')
                    a=statistic_class()
                    a.limma()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='LIMMA',message='Microarrays must be preprocessed before quality visualizations')
            def two_channel_segmentation(event):
                if workflow.count('tidy_up')>=1:
                    workflow.append('segmentation')
                    a=statistic_class()
                    a.segmentation_process()
                if workflow.count('tidy_up')==0:
                    warning=showinfo(title='Segmentation',message='Genome should be reordered before Segmentation')
            def two_channel_chrominfo(event):
                if workflow.count('preprocessing')>=1:
                    workflow.append('tidy_up')
                    analysis_type.append('tidy_up')
                    data_table()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Tidy Up',message='Microarrays must be preprocessed before reordeing')
            def two_channel_annotation(event):
                if workflow.count('limma')>=1:
                    a=functional_analysis()
                    a.retrieve_annotations()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Differential expression must be applied before gene annotation')
            def two_channel_clustering(event):
                if workflow.count('limma')>=1:
                    a=preprocess_class()
                    a.hierarchical_clustering()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
            if genomic_technology[len(genomic_technology)-1]=='Genepix':
                Label(win,text='GENEPIX DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            if genomic_technology[len(genomic_technology)-1]=='Agilent':
                Label(win,text='AGILENT DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                Label(win,text='AGILENT CGH DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            if genomic_technology[len(genomic_technology)-1]=='Imagene':
                Label(win,text='IMAGENE DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            if genomic_technology[len(genomic_technology)-1]=='Scanarray':
                Label(win,text='SCANARRAY DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            if genomic_technology[len(genomic_technology)-1]=='Quantarray':
                Label(win,text='QUANTARRAY DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            preprocess=Label(win,text='1. Microarrays Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            preprocess.place(x=7,y=100)
            preprocess.bind("<Button-1>",summarize_two_channel)
            quality=Label(win,text='2. Quality Control',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            quality.place(x=7,y=130)
            quality.bind("<Button-1>",two_channel_quality)
            normalization=Label(win,text='3. Normalization',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            normalization.place(x=7,y=160)
            normalization.bind("<Button-1>",two_channel_normalization)
            if genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                differential_expression=Label(win,text='4. Expression Analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
                differential_expression.place(x=7,y=190)
                differential_expression.bind("<Button-1>",two_channel_limma)
            if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                tidy_up=Label(win,text='4. Tidy Up',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
                tidy_up.place(x=7,y=190)
                tidy_up.bind("<Button-1>",two_channel_chrominfo)
            if genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                functional=Label(win,text='5. Annotation Retriever',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
                functional.place(x=7,y=220)
                functional.bind("<Button-1>", two_channel_annotation)
            if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                segmentation=Label(win,text='5. Segmentation',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
                segmentation.place(x=7,y=220)
                segmentation.bind("<Button-1>",two_channel_segmentation)
            if genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                Cluster=Label(win,text='6. Hierarchical Clustering',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
                Cluster.place(x=7,y=250)
                Cluster.bind("<Button-1>",two_channel_clustering)
            win.update()
        def qPCR(self):
            presentation_text.destroy()  
            def summarize_qPCR(event):
                if workflow.count('preprocessing')>0:
                    warning=showinfo(title='Workflow',message='Ct values already preprocessed.')
                if workflow.count('preprocessing')==0:
                    workflow.append('preprocessing')
                    data_table()
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- qPCR Ct values imported. Normalization is highly recommended.\n')
                    information_text.config(state=DISABLED)
            def qPCR_quality(event):
                if workflow.count('preprocessing')>=1:
                    a=quality_class()
                    a.array_quality()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Ct values must be preprocessed before quality visualizations')
            def qPCR_imputation(event):
                if workflow.count('preprocessing')>=1:
                    a=preprocess_class()
                    a.data_imputation()
                    workflow.append('imputation')
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Ct values must be preprocessed before quality visualizations')
            def qPCR_normalization(event):
                if workflow.count('preprocessing')>=1:
                    a=preprocess_class()
                    a.data_normalization()
                    workflow.append('normalization')
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Ct values must be preprocessed before quality visualizations')
            def qPCR_limma(event):
                if workflow.count('preprocessing')>=1:
                    workflow.append('limma')
                    a=statistic_class()
                    a.limma()
                if workflow.count('preprocessing')==0:
                    warning=showinfo(title='Workflow',message='Ct values must be preprocessed before quality visualizations')
            def qPCR_clustering(event):
                if workflow.count('limma')>=1:
                    a=preprocess_class()
                    a.hierarchical_clustering()
                if workflow.count('limma')==0:
                    warning=showinfo(title='Workflow',message='Hierarchical clustering needs from LiMMA or SAM analysis to be applied first.')
            Label(win,text='qPCR DATA',width=22,font=font_type2,bg='#8080c0',fg='#ffffff',relief='groove').place(x=7,y=37)
            preprocess=Label(win,text='1. Ct values Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            preprocess.place(x=7,y=100)
            preprocess.bind("<Button-1>",summarize_qPCR)
            quality=Label(win,text='2. Quality Control',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            quality.place(x=7,y=130)
            quality.bind("<Button-1>",qPCR_quality)
            imputation=Label(win,text='3. Data imputation',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            imputation.place(x=7,y=160)
            imputation.bind("<Button-1>",qPCR_imputation)
            normalization=Label(win,text='4. Normalization',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            normalization.place(x=7,y=190)
            normalization.bind("<Button-1>",qPCR_normalization)
            differential_expression=Label(win,text='5. Differential Expression',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            differential_expression.place(x=7,y=220)
            differential_expression.bind("<Button-1>",qPCR_limma)
            Cluster=Label(win,text='6. Hierarchical Clustering',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2',relief='flat')
            Cluster.place(x=7,y=250)
            Cluster.bind("<Button-1>", qPCR_clustering)
            win.update()
    #Funcion para abrir un nuevo proyecto
    def start_new_project():
        warning=askquestion(title='New Project',message='All data will be lost after creating a new project\nDo you want to continue?')
        if warning=='yes':
            new_project()
    #Funcion para generar tabla de datos
    def data_table():
        c=8;d=6
        if analysis_type[len(analysis_type)-1]=='none':
            for i in range(c):
                cols = []
                for j in range(d):
                    e=f[i]+'_'+f[j]
                    e=Entry(win,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                    e.place(x=320+x_resolution/2+j*109, y=45+20*i)
                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
        if analysis_type[len(analysis_type)-1]=='raw':
            for i in range(c):
                cols = []
                for j in range(d):
                    e=f[i]+'_'+f[j]
                    e=Entry(win,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                    e.place(x=320+x_resolution/2+j*109, y=45+20*i)
                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
            message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            for i in range(0,len(tables),1):
                tables[i].delete(0,END)
            if file_type[len(file_type)-1]=='txt':
                if data_transpose[len(data_transpose)-1]=='no':
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                if data_transpose[len(data_transpose)-1]=='yes':
                    infile=r.read_table('c:\\PG\\processes\\data\\custom_result.txt',sep='\t')
            if file_type[len(file_type)-1]=='vsc':
                if data_transpose[len(data_transpose)-1]=='no':
                    infile=r.read_table(file_name[len(file_name)-1],sep=',')
                if data_transpose[len(data_transpose)-1]=='yes':
                    infile=r.read_table('c:\\PG\\processes\\data\\custom_result.txt',sep='\t')
            infile=r.matrix(infile)
            b=infile
            d=len(b)
            z=len(samples)
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    if file_type[len(file_type)-1]=='txt':
                        infile=r.read_table('c:\\PG\\processes\\data\\custom_result.txt',sep='\t',col_names=col_name[0:d])
                    if file_type[len(file_type)-1]=='vsc':
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Data importing',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='annotations':
            for i in range(0,len(tables),1):
                tables[i].delete(0,END)
            message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            b=infile
            d=len(b)
            z=len(samples)
            if d<=8:
                for i in range(1,len(annotate_list)+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,len(annotate_list)+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Data importing',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='script':
            try:
                message_label=Label(win,text='Analyzing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                message_label.place(x=10,y=550)
                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                win.update()
                infile=open('c:\\PG\\processes\\data\\samples.txt','w')
                for i in range(0,len(samples),1):
                    infile.write(samples[i]+'\n')
                infile.close()
                r.source(script_name[len(script_name)-1])
                r('a<-'+function_name[len(function_name)-1]+'("c:/PG/processes/data/result.txt")')
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=len(samples)
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Custom Script',message='Impossible to apply function. No method selected')
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Custom Script',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='normalization':
            message_label=Label(win,text='Normalizing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            infile=open('c:\\PG\\processes\\data\\samples.txt','w')
            for i in range(0,len(samples),1):
                infile.write(samples[i]+'\n')
            infile.close()
            if column_normalization.get()==0:
                normalization_method='scale'
            if column_normalization.get()==1:
                normalization_method='quantile'
            if column_normalization.get()==2:
                normalization_method='z_score'
            if column_normalization.get()==3:
                normalization_method='fold_change'
            r.source("c:\\PG\\srs.R")
            try:
                if len(analysis_type)<=3:
                    if analysis_type[len(analysis_type)-2]=='Affymetrix':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-2]=='All-Affymetrix':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-2]=='Illumina':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-2]=='qPCR' or analysis_type[len(analysis_type)-2]=='imputate':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-2]=='raw':
                        if file_type[len(file_type)-1]=='txt':
                            r.normalization('c:\\PG\\processes\\data\\custom_result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                        if file_type[len(file_type)-1]=='vsc':
                            r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-2]=='Two channel' or analysis_type[len(analysis_type)-2]=='CGH' :
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                        r.anormalization('c:\\PG\\processes\\data\\box_result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if active_workflow[len(active_workflow)-1]!='none':
                        if genomic_technology[len(genomic_technology)-1]=='Custom':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                        if genomic_technology[len(genomic_technology)-1]!='Custom' and  genomic_technology[len(genomic_technology)-1]!='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    work_on_value.set(1)
                    if user_workflow_status[len(user_workflow_status)-1]!='open':
                        warning=showinfo(title='Normalization',message='Check Quality Control step again to see data corrections applied')
                if len(analysis_type)>3:
                    if analysis_type[len(analysis_type)-3]=='Affymetrix':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-3]=='All-Affymetrix':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-3]=='Illumina':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-3]=='qPCR' or analysis_type[len(analysis_type)-2]=='imputate':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-3]=='raw':
                        r.normalization('c:\\PG\\processes\\data\\custom_result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if analysis_type[len(analysis_type)-3]=='Two channel' or analysis_type[len(analysis_type)-3]=='CGH':
                        r.normalization('c:\\PG\\processes\\data\\result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                        r.anormalization('c:\\PG\\processes\\data\\box_result.txt',normalization_method,baseline_variable[len(baseline_variable)-1])
                    if active_workflow[len(active_workflow)-1]!='none':
                        if genomic_technology[len(genomic_technology)-1]=='Custom':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                        if genomic_technology[len(genomic_technology)-1]!='Custom' and  genomic_technology[len(genomic_technology)-1]!='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    work_on_value.set(1)
                    if user_workflow_status[len(user_workflow_status)-1]!='open':
                        warning=showinfo(title='Normalization',message='Check Quality Control step again to see data corrections applied')
                information_text.config(state=NORMAL)
                information_text.insert(INSERT,'NORMALIZATION\n')
                if column_normalization.get()==0:
                    information_text.insert(INSERT,'- Scale normalization method selected.\n')
                if column_normalization.get()==1:
                    information_text.insert(INSERT,'- Quantile normalization method selected.\n')
                if column_normalization.get()==2:
                    information_text.insert(INSERT,'- Z-Score normalization method selected.\n')
                if column_normalization.get()==3:
                    information_text.insert(INSERT,'- Fold change normalization method selected.\n')
                information_text.config(state=DISABLED)                        
            except:
                warning=showwarning(title='Normalization',message='Impossible to Normalize data')
                if active_workflow[len(active_workflow)-1]!='none':
                    if genomic_technology[len(genomic_technology)-1]!='Custom' and genomic_technology[len(genomic_technology)-1]!='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                    if genomic_technology[len(genomic_technology)-1]=='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=len(samples)
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Normalization',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='tidy_up':
            message_label=Label(win,text='Tidying Up...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            r.source("c:\\PG\\srs.R")
            def genome_view(event):
                r.tidyview('c:\\PG\\processes\\data\\result.txt')
            if workflow.count('tidy_up')>=1:
                try:
                    if workflow.count('normalization')>=1:
                        r.chrominfo('c:\\PG\\processes\\data\\result.txt')
                    if workflow.count('normalization')==0:
                        r.chrominfo('c:\\PG\\processes\\data\\raw_data.txt')
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    tidy_up_view=Label(win,text='view',font=font_type,bg='#fcfcfc',fg='#8080c0',relief='flat',cursor='hand2')
                    tidy_up_view.place(x=120,y=195)
                    tidy_up_view.bind("<Button-1>",genome_view)
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'CHOMOSOME INFORMATION\n')
                    information_text.insert(INSERT,'- Chromosomes and positions added.\n')
                    information_text.config(state=DISABLED)
                except:
                    warning=showwarning(title='Tidy Up',message='Impossible to reorder information')
                    Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=len(samples)
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                try:
                    if d>8:
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                        infile=r.matrix(infile)
                        b=infile
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                except:
                    warning=showwarning(title='Tidy Up',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='segmentation':
            message_label=Label(win,text='Segmentation...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            if workflow.count('tidy_up')>=1:
                if workflow.count('normalization')>=1:
                    try:
                        r.source("c:\\PG\\srs.R")
                        r.segmentation('c:\\PG\\processes\\data\\result.txt','c:\\PG\\processes\\data\\normalized_a_values.txt','c:\\PG\\processes\\data\\design.txt',segmentation_list[len(segmentation_list)-1])
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                        information_text.config(state=NORMAL)
                        information_text.insert(INSERT,'GENOME INFORMATION\n')
                        information_text.insert(INSERT,'- Genome information reordered, gain and loss possitions provided.\n')
                        information_text.config(state=DISABLED)
                    except:
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                if workflow.count('normalization')==0:
                    try:
                        r.source("c:\\PG\\srs.R")
                        r.segmentation('c:\\PG\\processes\\data\\result.txt','c:\\PG\\processes\\data\\raw_data.txt','c:\\PG\\processes\\data\\design.txt',segmentation_list[len(segmentation_list)-1])
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                        information_text.config(state=NORMAL)
                        information_text.insert(INSERT,'GENOME INFORMATION\n')
                        information_text.insert(INSERT,'- Genome information reordered, gain and loss possitions provided.\n')
                        information_text.config(state=DISABLED)
                    except:
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=len(samples)
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                try:
                    if d>8:
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                        infile=r.matrix(infile)
                        b=infile
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                except:
                    warning=showwarning(title='Segmentation',message='Number of columns exceeded limit supported in GeneQuery Basic version')
        if analysis_type[len(analysis_type)-1]=='imputate':
            try:
                message_label=Label(win,text='Imputing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                message_label.place(x=10,y=550)
                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                win.update()
                infile=open('c:\\PG\\processes\\data\\samples.txt','w')
                for i in range(0,len(samples),1):
                    infile.write(samples[i]+'\n')
                infile.close()
                r.source("c:\\PG\\srs.R")
                try:
                    if len(analysis_type)<=3:
                        if genomic_technology[len(genomic_technology)-1]!='qPCR':
                            if file_type[len(file_type)-1]=='txt':
                                r.imputation('c:\\PG\\processes\\data\\result.txt',replaces.get(),replace_values.get())
                            if file_type[len(file_type)-1]=='vsc':
                                r.imputation('c:\\PG\\processes\\data\\result.txt',replaces.get(),replace_values.get())
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            r.imputation('c:\\PG\\processes\\data\\result.txt',replaces.get(),replace_values.get())
                    if len(analysis_type)>3:
                        r.imputation('c:\\PG\\processes\\data\\result.txt',replaces.get(),replace_values.get())
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- Data imputed by '+ replaces.get()+'\n')
                    information_text.config(state=DISABLED)
                except:
                    warning=showwarning(title='Imputation',message='Impossible to apply imputation on data')
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=len(samples)
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                try:
                    if d>8:
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                        infile=r.matrix(infile)
                        b=infile
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                except:
                    warning=showwarning(title='Imputation',message='Number of columns exceeded limit supported in Desktop version')
                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
            except:
                warning=showwarning(title='Imputation',message='Impossible to apply imputation on data')
                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                    Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=165)
        if analysis_type[len(analysis_type)-1]=='summarize':
            message_label=Label(win,text='Summarizing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            infile=open('c:\\PG\\processes\\data\\samples.txt','w')
            for i in range(0,len(samples),1):
                infile.write(samples[i]+'\n')
            infile.close()
            r.source("c:\\PG\\srs.R")
            if len(analysis_type)<=3:
                if file_type[len(file_type)-1]=='txt':
                    try:
                        r.summarization(file_name[len(file_name)-1],summarization_replaces.get())
                        information_text.config(state=NORMAL)
                        information_text.insert(INSERT,'- Data summarized by '+summarization_replaces.get()+'\n')
                        information_text.config(state=DISABLED)
                    except:
                        warning=showwarning(title='Summarize',message='Impossible to summarize data. Please, check your dataset')
                if file_type[len(file_type)-1]=='vsc':
                    try:
                        r.summarization('c:\\PG\\processes\\data\\result.txt',summarization_replaces.get())
                        information_text.config(state=NORMAL)
                        information_text.insert(INSERT,'- Data summarized by '+summarization_replaces.get()+'\n')
                        information_text.config(state=DISABLED)
                    except:
                        warning=showwarning(title='Summarize',message='Impossible to summarize data. Please, check the dataset')
            if len(analysis_type)>3:
                try:
                    r.summarization('c:\\PG\\processes\\data\\result.txt',summarization_replaces.get())
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- Data summarized by '+summarization_replaces.get()+'\n')
                    information_text.config(state=DISABLED)
                except:
                    warning=showwarning(title='Summarize',message='Impossible to summarize data. Please, check the dataset')
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=len(samples)
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Summarize',message='Number of columns exceeded limit supported in Desktop version')
        if analysis_type[len(analysis_type)-1]=='limma':
            message_label=Label(win,text='LIMMA analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            infile=open('c:\\PG\\processes\\data\\samples.txt','w')
            for i in range(0,len(samples),1):
                infile.write(samples[i]+'\n')
            infile.close()
            r.source("c:\\PG\\srs.R")
            if fdr.get()=='Benjamini-Hochberg':
                method='BH'
            if fdr.get()=='Holm':
                method='holm'
            if fdr.get()=='Bonferroni':
                method='bonferroni'
            if fdr.get()=='none':
                method='none'
            if len(analysis_type)>3:
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    r.limma('All-Affymetrix','c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='Illumina':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                workflow.append('limma')
            if len(analysis_type)<=3:
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    r.limma('All-Affymetrix','c:\\PG\\processes\\data\\result.txt',method)
                if analysis_type[len(analysis_type)-2]=='Illumina' :
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)   
                if analysis_type[len(analysis_type)-2]=='raw' or genomic_technology[len(genomic_technology)-1]=='Custom':
                    if file_type[len(file_type)-1]=='txt':
                        r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                    if file_type[len(file_type)-1]=='vsc':
                        r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    r.limma(genomic_technology[len(genomic_technology)-1],'c:\\PG\\processes\\data\\result.txt',method)
                workflow.append('limma')
            information_text.config(state=NORMAL)
            information_text.insert(INSERT,'LINEAR MODELS FOR MICROARRAYS.\n')
            if fdr.get()=='Benjamini-Hochberg':
                information_text.insert(INSERT,'- Benjamini-Hochberg method selected for False Discovery Rate correction.\n')
            if fdr.get()=='Holm':
                information_text.insert(INSERT,'- Holm method selected for False Discovery Rate correction.\n')
            if fdr.get()=='Bonferroni':
                information_text.insert(INSERT,'- Bonferroni method selected for False Discovery Rate correction.\n')
            if fdr.get()=='none':
                information_text.insert(INSERT,'- No False Discovery Rate (FDR) correction applied. FDR is highly recommended to get confident results.\n')
            information_text.config(state=DISABLED)
            if active_workflow[len(active_workflow)-1]!='none':
                if genomic_technology[len(genomic_technology)-1]!='Custom' and genomic_technology[len(genomic_technology)-1]!='qPCR':  
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                if genomic_technology[len(genomic_technology)-1]=='Custom':
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=6
            if d<=8:
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                if genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                    if columns[0]=='GeneID' or columns[0]!='ProbeName'or columns[0]=='GeneName':
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                    if columns[0]!='GeneID' and columns[0]=='ProbeName':
                        for i in range(1,z+2,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    for i in range(1,z+2,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                    if genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                        if columns[0]=='GeneID' or columns[0]=='ProbeName' or columns[0]=='GeneName':
                            for i in range(1,z+1,1):
                                if columns.count(b[d-i][0][0])==0:
                                    columns.append(b[d-i][0][0])
                        if columns[0]!='GeneID' and columns[0]!='ProbeName':
                            for i in range(1,z+2,1):
                                if columns.count(b[d-i][0][0])==0:
                                    columns.append(b[d-i][0][0])
                    if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        for i in range(1,z+2,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='LIMMA',message='Number of columns exceeded limit supported in Desktop version')
                if active_workflow[len(active_workflow)-1]!='none':
                    if genomic_technology[len(genomic_technology)-1]!='Custom' and genomic_technology[len(genomic_technology)]!='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if genomic_technology[len(genomic_technology)-1]=='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
        if analysis_type[len(analysis_type)-1]=='ttest':
            message_label=Label(win,text='T-TEST analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            infile=open('c:\\PG\\processes\\data\\samples.txt','w')
            for i in range(0,len(samples),1):
                infile.write(samples[i]+'\n')
            infile.close()
            r.source("c:\\PG\\srs.R")
            if fdr.get()=='Benjamini-Hochberg':
                method='BH'
            if fdr.get()=='Holm':
                method='holm'
            if fdr.get()=='Bonferroni':
                method='bonferroni'
            if fdr.get()=='none':
                method='none'
            try:
                if len(analysis_type)>3:
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    workflow.append('limma')
                if len(analysis_type)<=3:
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if analysis_type[len(analysis_type)-2]=='Illumina':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')   
                    if analysis_type[len(analysis_type)-2]=='raw' or genomic_technology[len(genomic_technology)-1]=='Custom':
                        if file_type[len(file_type)-1]=='txt':
                            r.ttest('c:\\PG\\processes\\data\\result.txt')
                        if file_type[len(file_type)-1]=='vsc':
                            r.ttest('c:\\PG\\processes\\data\\result.txt')
                    if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        r.ttest('c:\\PG\\processes\\data\\result.txt')
                    workflow.append('limma')
                information_text.config(state=NORMAL)
                information_text.insert(INSERT,'UNPAIRED T-TEST.\n')
                if fdr.get()=='Benjamini-Hochberg':
                    information_text.insert(INSERT,'- Benjamini-Hochberg method selected for False Discovery Rate correction.\n')
                if fdr.get()=='Holm':
                    information_text.insert(INSERT,'- Holm method selected for False Discovery Rate correction.\n')
                if fdr.get()=='Bonferroni':
                    information_text.insert(INSERT,'- Bonferroni method selected for False Discovery Rate correction.\n')
                if fdr.get()=='none':
                    information_text.insert(INSERT,'- No False Discovery Rate (FDR) correction applied. FDR is highly recommended to get confident results.\n')
                information_text.config(state=DISABLED)
            except:
                warning=showwarning(title='SAM',message='Impossible to apply T-Test analysis')
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=4
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                    if genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Unpaired T-Test',message='Number of columns exceeded limit supported in Desktop version')
                if active_workflow[len(active_workflow)-1]!='none':
                    if genomic_technology[len(genomic_technology)-1]!='Custom' and genomic_technology[len(genomic_technology)-1]!='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if genomic_technology[len(genomic_technology)-1]=='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
        if analysis_type[len(analysis_type)-1]=='sam':
            message_label=Label(win,text='SAM analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)                
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            infile=open('c:\\PG\\processes\\data\\samples.txt','w')
            for i in range(0,len(samples),1):
                infile.write(samples[i]+'\n')
            infile.close()
            r.source("c:\\PG\\srs.R")
            if statistic.get()=='Standard':
                statistics='standard'
            if statistic.get()=='Wilcoxon':
                statistics='wilcoxon'
            if fdr.get()=='Benjamini-Hochberg':
                method='BH'
            if fdr.get()=='Holm':
                method='holm'
            if fdr.get()=='Bonferroni':
                method='bonferroni'
            if fdr.get()=='none':
                method='none'
            try:
                if len(analysis_type)>3:
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        r.sam('c:\\PG\\processes\\data\\result.txt','Affymetrix',method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='qPCR':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    workflow.append('limma')
                if len(analysis_type)<=3:
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                        r.sam('c:\\PG\\processes\\data\\result.txt','Affymetrix',method,statistics)
                    if analysis_type[len(analysis_type)-2]=='Illumina':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if analysis_type[len(analysis_type)-2]=='raw' or genomic_technology[len(genomic_technology)-1]=='Custom':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        r.sam('c:\\PG\\processes\\data\\result.txt',genomic_technology[len(genomic_technology)-1],method,statistics)
                    workflow.append('limma')
                information_text.config(state=NORMAL)
                information_text.insert(INSERT,'SIGNIFICANT ANALYSIS OF MICROARRAYS.\n')
                if statistic.get()=='Standard':
                    information_text.insert(INSERT,'- Statistic method selected: Standard.\n')
                if statistic.get()=='Wilcoxon':
                    information_text.insert(INSERT,'- Statistic method selected: Wilcoxon.\n')
                if fdr.get()=='Benjamini-Hochberg':
                    information_text.insert(INSERT,'- Benjamini-Hochberg method selected for False Discovery Rate correction.\n')
                if fdr.get()=='Holm':
                    information_text.insert(INSERT,'- Holm method selected for False Discovery Rate correction.\n')
                if fdr.get()=='Bonferroni':
                    information_text.insert(INSERT,'- Bonferroni method selected for False Discovery Rate correction.\n')
                if fdr.get()=='none':
                    information_text.insert(INSERT,'- No False Discovery Rate (FDR) correction applied. FDR is highly recommended to get confident results.\n')
                information_text.config(state=DISABLED)
            except:
                warning=showwarning(title='SAM',message='Impossible to apply SAM analysis')
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=6
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='SAM',message='Number of columns exceeded limit supported in Desktop version')
                if active_workflow[len(active_workflow)-1]!='none':
                    if genomic_technology[len(genomic_technology)-1]!='Custom' and genomic_technology[len(genomic_technology)-1]!='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if genomic_technology[len(genomic_technology)-1]=='qPCR':
                        Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=225)
            process.create_line(0,7,200,7,fill='#8080c0',width=15)
            win.update()
        if analysis_type[len(analysis_type)-1]=='Affymetrix' or analysis_type[len(analysis_type)-1]=='All-Affymetrix' :
            message_label=Label(win,text='AFFYMETRIX analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            r.source('c:\\PG\\srs.R')
            if analysis_type[len(analysis_type)-1]=='Affymetrix' or analysis_type[len(analysis_type)-1]=='All-Affymetrix':
                try:
                    if affymetrix_method[len(affymetrix_method)-1]=='MAS5':
                        r.affymetrix('c:\\PG\\processes\\data\\Affymetrix\\design.txt',affymetrix_method[len(affymetrix_method)-1])
                    if affymetrix_method[len(affymetrix_method)-1]=='RMA':
                        r.affymetrix('c:\\PG\\processes\\data\\Affymetrix\\design.txt',affymetrix_method[len(affymetrix_method)-1])                        
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                    process.create_line(0,7,100,7,fill='#8080c0',width=15)
                    win.update()
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                    infile=r.matrix(infile)
                    process.create_line(0,7,150,7,fill='#8080c0',width=15)
                    win.update()
                    b=infile
                    d=len(b)
                    z=6
                    if d<=8:
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                    try:
                        if d>8:
                            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                            infile=r.matrix(infile)
                            b=infile
                            for i in range(1,z+1,1):
                                if columns.count(b[d-i][0][0])==0:
                                    columns.append(b[d-i][0][0])
                    except:
                        warning=showwarning(title='Affymetrix',message='Number of columns exceeded limit supported in Desktop version')
                    process.create_line(0,7,200,7,fill='#8080c0',width=15)
                    win.update()
                except:
                    Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                    warning=showwarning(title='Affymetrix',message='Affymetrix chip not supported by GeneQuery Basic edition')
                    process.delete(ALL)
                    message_label.destroy()
                    new_project()
        if analysis_type[len(analysis_type)-1]=='Two channel' or analysis_type[len(analysis_type)-1]=='CGH' :
            if genomic_technology[len(genomic_technology)-1]=='Genepix':
                message_label=Label(win,text='GENEPIX analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                technology='genepix'
            if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' :
                message_label=Label(win,text='AGILENT analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                technology='agilent'
            if genomic_technology[len(genomic_technology)-1]=='Imagene':
                message_label=Label(win,text='IMAGENE analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                technology='imagene'
            if genomic_technology[len(genomic_technology)-1]=='Scanarray':
                message_label=Label(win,text='SCANARRAY analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                technology='scanarray'
            if genomic_technology[len(genomic_technology)-1]=='Quantarray':
                message_label=Label(win,text='QUANTARRAY analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
                technology='quantarray'
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            r.source('c:\\PG\\srs.R')
            if genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                r.twochannel('c:\\PG\\processes\\data\\twochannel\\design.txt',technology,background.get(),within_normalization.get(),'')
            if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                r.twochannel('c:\\PG\\processes\\data\\twochannel\\design.txt',technology,background.get(),within_normalization.get(),'CGH')
            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
            process.create_line(0,7,100,7,fill='#8080c0',width=15)
            win.update()
            infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
            infile=r.matrix(infile)
            process.create_line(0,7,150,7,fill='#8080c0',width=15)
            win.update()
            b=infile
            d=len(b)
            z=6
            if d<=8:
                for i in range(1,z+1,1):
                    if columns.count(b[d-i][0][0])==0:
                        columns.append(b[d-i][0][0])
            try:
                if d>8:
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                    infile=r.matrix(infile)
                    b=infile
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
            except:
                warning=showwarning(title='Two Channel',message='Number of columns exceeded limit supported in Desktop version')
            process.create_line(0,7,200,7,fill='#8080c0',width=15)
            win.update()
        if analysis_type[len(analysis_type)-1]=='Illumina':
            skip=[0]
            message_label=Label(win,text='ILLUMINA analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            file=open('c:\\PG\\processes\\data\\Illumina\\Design.txt','r')
            content=file.readline()
            n1=content.find('\t')
            file_names=content[n1+1:len(content)-1]
            files=open(file_names,'r')
            contents=files.readlines()
            for i in range(0,len(contents),1):
                if contents[i].count('ProbeID')>0 and contents[i].count('AVG_Signal')>0:
                    skip.append(i)
            r.source('c:\\PG\\srs.R')
            try:
                r.illumina(file_names,skip[len(skip)-1])
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=6
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                try:
                    if d>8:
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                        infile=r.matrix(infile)
                        b=infile
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                except:
                    warning=showwarning(title='Illumina',message='Number of columns exceeded limit supported in Basic version')
                Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                process.create_line(0,7,200,7,fill='#8080c0',width=15)
                win.update()
            except:
                warning=showwarning(title='Illumina',message='qPCR data format not supported by GeneQuery Basic edition')
                Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                new_project()
        if analysis_type[len(analysis_type)-1]=='qPCR':
            message_label=Label(win,text='qPCR analysis running...',font=font_type3,bg='#fcfcfc',relief='flat')
            message_label.place(x=10,y=550)
            process.create_line(0,7,50,7,fill='#8080c0',width=15)
            win.update()
            file=open('c:\\PG\\processes\\data\\qPCR\\Design.txt','r')
            content=file.readline()
            n1=content.find('\t')
            file_names=content[n1+1:len(content)-1]
            files=open(file_names,'r')
            contents=files.readlines()
            r.source('c:\\PG\\srs.R')
            try:
                r.qPCR(file_names,aggregate_method[len(aggregate_method)-1],qpcr_missing[len(qpcr_missing)-1])
                process.create_line(0,7,100,7,fill='#8080c0',width=15)
                win.update()
                infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                infile=r.matrix(infile)
                process.create_line(0,7,150,7,fill='#8080c0',width=15)
                win.update()
                b=infile
                d=len(b)
                z=6
                if d<=8:
                    for i in range(1,z+1,1):
                        if columns.count(b[d-i][0][0])==0:
                            columns.append(b[d-i][0][0])
                try:
                    if d>8:
                        infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t',col_names=col_name[0:d])
                        infile=r.matrix(infile)
                        b=infile
                        for i in range(1,z+1,1):
                            if columns.count(b[d-i][0][0])==0:
                                columns.append(b[d-i][0][0])
                except:
                    warning=showwarning(title='qPCR',message='Number of columns exceeded limit supported in Basic version')
                Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                process.create_line(0,7,200,7,fill='#8080c0',width=15)
                win.update()
            except:
                warning=showwarning(title='qPCR',message='qPCR chip not supported by GeneQuery Basic edition')
                Label(win,image=foto18,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                new_project()
        if analysis_type[len(analysis_type)-1]!='none'and analysis_type[len(analysis_type)-1]!='anova':
            if  c >8:
                c=8
            if d > 6:
                d=6
            rows=[]
            for i in range(c):
                try:
                    cols = []
                    for j in range(d):
                        e=Entry(win,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                        e.place(x=320+x_resolution/2+j*109, y=45+20*i)
                        e.config(state=NORMAL,disabledbackground='#fcfcfc')
                        e.insert(INSERT,b[j][0][i])
                        e.config(state=DISABLED,disabledbackground='#fcfcfc')
                        cols.append(e)
                except:
                    if analysis_type[len(analysis_type)-1]!='raw':
                        e.config(state=NORMAL,disabledbackground='#fcfcfc')
                        e.insert(INSERT,'GeneID')
                        e.config(state=DISABLED,disabledbackground='#fcfcfc')
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                    infile=r.matrix(infile)
                    break
                rows.append(cols)
            message_label.destroy()
            process.create_line(0,7,200,7,fill='#8080c0',width=15)
            win.update()
            process.delete(ALL)
            def scroll_move(event):
                def move_table():
                    b=infile
                    c=8
                    d=6
                    rows=[]
                    for i in range(c):
                        cols = []
                        for j in range(d):
                            e=Entry(win,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                            e.place(x=320+x_resolution/2+j*109, y=45+20*i)
                            if event.x<100:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.x>100 and event.x<200:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+6][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.x>200 and event.x<300:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+12][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.x>300 and event.x<400:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+18][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.x>400 and event.x<500:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+24][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.x>500 and event.x<600:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+30][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            cols.append(e)
                        rows.append(cols)
                scroll_table.delete(ALL)
                scroll_position.append(event.x)
                if event.x>=0 and event.x<=655:
                    scrolling_bar=scroll_table.create_rectangle(event.x,3,event.x+30,5,fill='#8080c0')
                    move_table()
                if event.x<0:
                    scrolling_bar=scroll_table.create_rectangle(3,3,30,5,fill='#8080c0')
                if event.x>628:
                    scrolling_bar=scroll_table.create_rectangle(628,3,655,5,fill='#8080c0')
                win.update()
            def scroll_move2(event):
                def move_table():
                    b=infile
                    c=8
                    d=6
                    rows=[]
                    if scroll_position[len(scroll_position)-1]<=100:
                        x=0
                    if scroll_position[len(scroll_position)-1]>100 and scroll_position[len(scroll_position)-1]<=200 :
                        x=6
                    if scroll_position[len(scroll_position)-1]>200 and scroll_position[len(scroll_position)-1]<=300 :
                        x=12
                    if scroll_position[len(scroll_position)-1]>300 and scroll_position[len(scroll_position)-1]<=400 :
                        x=18
                    if scroll_position[len(scroll_position)-1]>400 and scroll_position[len(scroll_position)-1]<=500 :
                        x=24
                    if scroll_position[len(scroll_position)-1]>500 and scroll_position[len(scroll_position)-1]<=600 :
                        x=30
                    for i in range(c):
                        cols = []
                        for j in range(d):
                            e=Entry(win,font=font_type,width=18,relief='groove',bg='#fcfcfc',fg='#8080c0')
                            e.place(x=320+x_resolution/2+j*109, y=45+20*i)
                            if event.y<20:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.y>20 and event.y<50:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i+8])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.y>50 and event.y<70:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i+16])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.y>70 and event.y<100:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i+24])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.y>100 and event.y<130:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i+32])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            if event.y>130 and event.y<160:
                                try:
                                    e.config(state=NORMAL,disabledbackground='#fcfcfc')
                                    e.insert(INSERT,b[j+x][0][i+40])
                                    e.config(state=DISABLED,disabledbackground='#fcfcfc')
                                except:
                                    pass
                            cols.append(e)
                        rows.append(cols)
                scroll_table2.delete(ALL)
                if event.y>=0 and event.y<=160:
                    scrolling_bar2=scroll_table2.create_rectangle(3,event.y,5,event.y+30,fill='#8080c0')
                    move_table()
                if event.y<0:
                    scrolling_bar2=scroll_table2.create_rectangle(3,3,5,30,fill='#8080c0')
                if event.y>160:
                    scrolling_bar2=scroll_table2.create_rectangle(3,133,5,160,fill='#8080c0')
                win.update()
            scroll_table=Canvas(win,width=655,height=5,bg='#fcfcfc',relief='groove',cursor='hand2')
            scroll_table.place(x=320+x_resolution/2,y=210)
            scroll_table2=Canvas(win,width=5,height=160,bg='#fcfcfc',relief='groove',cursor='hand2')
            scroll_table2.place(x=310+x_resolution/2,y=45)
            scrolling_bar=scroll_table.create_rectangle(3,3,30,5,fill='#8080c0')
            scrolling_bar2=scroll_table2.create_rectangle(3,3,5,30,fill='#8080c0')
            scroll_table.bind("<B1-Motion>",scroll_move)
            scroll_table2.bind("<B1-Motion>",scroll_move2)
            if analysis_type[len(analysis_type)-1]=='Affymetrix' or analysis_type[len(analysis_type)-1]=='All-Affymetrix' or analysis_type[len(analysis_type)-1]=='Illumina' or analysis_type[len(analysis_type)-1]=='Two channel' or analysis_type[len(analysis_type)-1]=='CGH' or analysis_type[len(analysis_type)-1]=='qPCR':
                del(columns[0:len(columns)])
                for i in range(0,len(infile)):
                    columns.append(infile[i][0][0])
                for i in range(1,len(infile)):
                    control_columns.append(infile[i][0][0])
                if analysis_type[len(analysis_type)-1]=='Affymetrix':
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                if analysis_type[len(analysis_type)-1]=='All-Affymetrix':
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(int(2**float(infile[1][0][i])));data2.append(int(2**float(infile[2][0][i])))
                if analysis_type[len(analysis_type)-1]=='Illumina':
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                if analysis_type[len(analysis_type)-1]=='Two channel' or analysis_type[len(analysis_type)-1]=='qPCR' or analysis_type[len(analysis_type)-1]=='CGH':
                    for i in range(1,len(infile[1][0]),1):
                        if infile[1][0][i]!='NA' and infile[2][0][i]!='NA':
                            gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                for i in range(0,len(data1),1):
                    if data1[i]>values[len(values)-1]:
                        values.append(data1[i])
                    if data2[i]>values2[len(values2)-1]:
                        values2.append(data2[i])
                if values2[len(values2)-1]>values[len(values)-1]:
                    values.append(values2[len(values2)-1])
                for i in range(0,len(data1),1):
                    new_data1.append((data1[i]*650)/values[len(values)-1])
                    new_data2.append((data2[i]*400)/values[len(values)-1])
                sample1.insert(INSERT,columns[2]);sample2.insert(INSERT,columns[1])
                graphic()
    #funcion analisis regresion lineal
    def regression(new_data1,new_data2):
        try:
            if len(new_data1) != len(new_data2):  raise ValueError, 'unequal length'
            N = len(new_data1)
            Sx = Sy = Sxx = Syy = Sxy = 0.0
            for x, y in map(None, new_data1, new_data2):
                Sx = Sx + x
                Sy = Sy + y
                Sxx = Sxx + x*x
                Syy = Syy + y*y
                Sxy = Sxy + x*y
            det = Sxx * N - Sx * Sx
            a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
            meanerror = residual = 0.0
            for x, y in map(None, new_data1, new_data2):
                meanerror = meanerror + (y - Sy/N)**2
                residual = residual + (y - a * x - b)**2
            RR = 1 - residual/meanerror
            ss = residual / (N-2)
            Var_a, Var_b = ss * N / det, ss * Sxx / det
            return a, b, RR
        except:
            showwarning(title='Regression',message='Impossible to apply regression test.')
            a='none'
            b='none'
            RR='none'
    #funcion para abrir archivos custom
    def graphic():
        graph_status=Canvas(win,width=100,height=10,bg='#fcfcfc',borderwidth=1,relief='flat')
        graph_status.place(x=580+x_resolution/2,y=430)
        graph_status.create_line(0,6,25,6,fill='#8080c0',width=10)
        win.update()
        graph_panel.delete(ALL)
        graph_panel.create_line(0,40,650,40,fill='#c5c5c5');graph_panel.create_line(0,80,650,80,fill='#c5c5c5')
        graph_panel.create_line(0,120,650,120,fill='#c5c5c5');graph_panel.create_line(0,160,650,160,fill='#c5c5c5')
        graph_panel.create_line(0,200,650,200,fill='#c5c5c5');graph_panel.create_line(0,240,650,240,fill='#c5c5c5')
        graph_panel.create_line(0,280,650,280,fill='#c5c5c5');graph_panel.create_line(0,320,650,320,fill='#c5c5c5')
        graph_panel.create_line(0,360,650,360,fill='#c5c5c5')
        graph_panel.create_line(65,0,65,400,fill='#c5c5c5');graph_panel.create_line(130,0,130,400,fill='#c5c5c5')
        graph_panel.create_line(195,0,195,400,fill='#c5c5c5');graph_panel.create_line(260,0,260,400,fill='#c5c5c5')
        graph_panel.create_line(325,0,325,400,fill='#c5c5c5');graph_panel.create_line(390,0,390,400,fill='#c5c5c5')
        graph_panel.create_line(455,0,455,400,fill='#c5c5c5');graph_panel.create_line(520,0,520,400,fill='#c5c5c5')
        graph_panel.create_line(585,0,585,400,fill='#c5c5c5')
        graph_status.create_line(0,6,50,6,fill='#8080c0',width=10)
        win.update()
        for i in range(0,len(new_data1),1):
            if shape[len(shape)-1]=='oval':
                a=graph_panel.create_oval(0,0,int(size[len(size)-1]),int(size[len(size)-1]),fill=color[len(color)-1])
            if shape[len(shape)-1]=='rectangle':
                a=graph_panel.create_rectangle(0,0,int(size[len(size)-1]),int(size[len(size)-1]),fill=color[len(color)-1])            
            if graph_type[len(graph_type)-1]=='linear':
                if graph_properties.get()=='Graph Options' or graph_properties.get()=='Zoom out' or graph_properties.get()=='Graph Settings' or graph_properties.get()=='Show regression values':
                    graph_panel.move(a,new_data1[i],400-(new_data2[i]))
                    if label_display.get()==1:
                        id=graph_panel.create_text(new_data1[i],390-(new_data2[i]),text=gene_id[i],font=(graph_font[len(graph_font)-1],graph_font_size[len(graph_font_size)-1]))
                if graph_properties.get()=='Zoom in':
                    graph_panel.move(a,3*(new_data1[i]),390-(3*(new_data2[i])))
                    if label_display.get()==1:
                        id=graph_panel.create_text(3*(new_data1[i]),390-(3*(new_data2[i])),text=gene_id[i],font=(graph_font[len(graph_font)-1],graph_font_size[len(graph_font_size)-1]))
        graph_status.create_line(0,6,75,6,fill='#8080c0',width=10)
        win.update()
        if genomic_technology[len(genomic_technology)-1]=='Custom' or genomic_technology[len(genomic_technology)-1]=='qPCR' or genomic_technology[len(genomic_technology)-1]=='Illumina' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
            if log_scale[len(log_scale)-1]=='TRUE':
                value1=float(values[len(values)-1]);value2=float(values[len(values)-2])
                if graph_properties.get()=='Zoom in':
                    value1=values[len(values)-1]/3;value2=values[len(values)-2]/3
            if log_scale[len(log_scale)-1]=='FALSE':
                value1=math.log(values[len(values)-1]);value2=math.log(values[len(values)-2])
                if graph_properties.get()=='Zoom in':
                    value1=math.log(values[len(values)-1])/3;value2=math.log(values[len(values)-2])/3
        if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
            if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                value1=values[len(values)-1];value2=values[len(values)-2]
                if graph_properties.get()=='Zoom in':
                    value1=values[len(values)-1]/3;value2=values[len(values)-2]/3
            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                value1=math.log(values[len(values)-1]);value2=math.log(values[len(values)-2])
                if graph_properties.get()=='Zoom in':
                    value1=math.log(values[len(values)-1]/3);value2=math.log(values[len(values)-2]/3)
        x1.config(state=NORMAL);x2.config(state=NORMAL);x3.config(state=NORMAL)
        x4.config(state=NORMAL);x5.config(state=NORMAL);x6.config(state=NORMAL)
        x7.config(state=NORMAL);x8.config(state=NORMAL);x9.config(state=NORMAL);x10.config(state=NORMAL)
        y1.config(state=NORMAL);y2.config(state=NORMAL);y3.config(state=NORMAL)
        y4.config(state=NORMAL);y5.config(state=NORMAL);y6.config(state=NORMAL)
        y7.config(state=NORMAL);y8.config(state=NORMAL);y9.config(state=NORMAL);y0.config(state=NORMAL)
        x10.delete(0,END);x10.insert(INSERT,value1)
        x9.delete(0,END);x9.insert(INSERT,value1*0.9);y9.delete(0,END);y9.insert(INSERT,value2*0.1)
        x8.delete(0,END);x8.insert(INSERT,value1*0.8);y8.delete(0,END);y8.insert(INSERT,value2*0.2)
        x7.delete(0,END);x7.insert(INSERT,value1*0.7);y7.delete(0,END);y7.insert(INSERT,value2*0.3)
        x6.delete(0,END);x6.insert(INSERT,value1*0.6);y6.delete(0,END);y6.insert(INSERT,value2*0.4)
        x5.delete(0,END);x5.insert(INSERT,value1*0.5);y5.delete(0,END);y5.insert(INSERT,value2*0.5)
        x4.delete(0,END);x4.insert(INSERT,value1*0.4);y4.delete(0,END);y4.insert(INSERT,value2*0.6)
        x3.delete(0,END);x3.insert(INSERT,value1*0.3);y3.delete(0,END);y3.insert(INSERT,value2*0.7)
        x2.delete(0,END);x2.insert(INSERT,value1*0.2);y2.delete(0,END);y2.insert(INSERT,value2*0.8)
        x1.delete(0,END);x1.insert(INSERT,value1*0.1);y1.delete(0,END);y1.insert(INSERT,value2*0.9)
        y0.delete(0,END);y0.insert(INSERT,value2)
        x1.config(state=DISABLED,disabledbackground='#ffffff');x2.config(state=DISABLED,disabledbackground='#ffffff');x3.config(state=DISABLED,disabledbackground='#ffffff')
        x4.config(state=DISABLED,disabledbackground='#ffffff');x5.config(state=DISABLED,disabledbackground='#ffffff');x6.config(state=DISABLED,disabledbackground='#ffffff')
        x7.config(state=DISABLED,disabledbackground='#ffffff');x8.config(state=DISABLED,disabledbackground='#ffffff');x9.config(state=DISABLED,disabledbackground='#ffffff');x10.config(state=DISABLED,disabledbackground='#ffffff')
        y1.config(state=DISABLED,disabledbackground='#ffffff');y2.config(state=DISABLED,disabledbackground='#ffffff');y3.config(state=DISABLED,disabledbackground='#ffffff')
        y4.config(state=DISABLED,disabledbackground='#ffffff');y5.config(state=DISABLED,disabledbackground='#ffffff');y6.config(state=DISABLED,disabledbackground='#ffffff')
        y7.config(state=DISABLED,disabledbackground='#ffffff');y8.config(state=DISABLED,disabledbackground='#ffffff');y9.config(state=DISABLED,disabledbackground='#ffffff');y0.config(state=DISABLED,disabledbackground='#ffffff')
        graph_status.create_line(0,6,100,6,fill='#8080c0',width=10)
        win.update()
        graph_status.destroy()
        regression_result=regression(new_data1,new_data2)
        slope.append(regression_result[0])
        intercept.append(regression_result[1])
        cor_value.append(regression_result[2])
        Y=regression_result[1]
        X=(-regression_result[1])/regression_result[0]
        graph_panel.create_line(0,400,650,0,fill='#ff0000',width=2)
    def export_main_graph():
        graph_panel.postscript(file='C:\\PG\\processes\\data\\main_graph.ps',colormode='color')
        try:
            os.system('start AcroRd32 C:\\PG\\processes\\data\\main_graph.ps')
        except:
            pass
    def show_regression_values():
        def quit_regression_values():
            regression_frame.destroy()
            label1.destroy();label_correlation.destroy()
            label2.destroy();label_adjusted.destroy()
            label3.destroy();correlation_entry.destroy()
            close.destroy()
            regression_window.append('close')
            graph_properties.set('Graph Options')
        if graph_properties.get()=='Show regression values':
            regression_frame=Frame(win,width=200,height=150,bg='#f5f5f5',relief='groove',borderwidth=2)
            regression_frame.place(x=325+x_resolution/2,y=245)
            label1=Label(win,text='     ',font=font_type,bg=color[len(color)-1],relief='groove')
            label1.place(x=335+x_resolution/2,y=255)
            label_adjusted=Label(win,text='Spot data', font=font_type,bg='#f5f5f5')
            label_adjusted.place(x=370+x_resolution/2,y=255)
            label2=Label(win,text='     ',font=font_type,bg='#ff0000',relief='groove')
            label2.place(x=335+x_resolution/2,y=285)
            label_correlation=Label(win,text='Perfect data line fit', font=font_type,bg='#f5f5f5')
            label_correlation.place(x=370+x_resolution/2,y=285)
            label3=Label(win,text='Correlation:',font=font_type,bg='#f5f5f5')
            label3.place(x=335+x_resolution/2,y=315)
            correlation_entry=Entry(win,width=10,font=font_type,bg='#f5f5f5',relief='flat')
            correlation_entry.place(x=400+x_resolution/2,y=315)
            close=Button(win,text='Close',width=7,bg='#ffffff',activebackground='#f5f5f5',font=font_type,relief='groove',cursor='hand2',command=quit_regression_values)
            close.place(x=400+x_resolution/2,y=350)
            regression_window.append('open')
            correlation_entry.delete(0,END)
            correlation_entry.insert(INSERT,str(cor_value[len(cor_value)-1]))
            correlation_entry.config(state=DISABLED,disabledbackground='#f5f5f5')
    def graphic_properties(event):
        if analysis_type[len(analysis_type)-1]!='none':
            if up_window_status[len(up_window_status)-1]=='close':
                def properties_selection(event):
                    if graph_properties.get()=='Graph Settings':
                        properties_menu.destroy()
                        graphic_options(event)
                        up_window_status.append('close')
                    if graph_properties.get()=='Zoom in':
                        properties_menu.destroy()
                        win.update()
                        up_window_status.append('close')
                        graphic()
                    if graph_properties.get()=='Zoom out':
                        properties_menu.destroy()
                        win.update()
                        up_window_status.append('close')
                        graphic()
                    if graph_properties.get()=='Export Current Graph':
                        properties_menu.destroy()
                        win.update()
                        up_window_status.append('close')
                        export_main_graph()
                    if graph_properties.get()=='Show regression values':
                        properties_menu.destroy()
                        win.update()
                        up_window_status.append('close')
                        show_regression_values()
                properties_menu=OptionMenu(win, graph_properties,'Zoom in', 'Zoom out','Show regression values','Graph Settings')
                properties_menu.place(x=position_x[len(position_x)-1]+320,y=position_y[len(position_y)-1]+240)
                properties_menu.bind("<Button-1>",properties_selection)
                graph_properties.set('Graph Options')
                up_window_status.append('open')
    def graphic_options(event):
        def cancel_graphic_options():
            top.destroy()
            window_status.append('close')
        def accept_options():
            color.append(color_types.get())
            shape.append(shape_types.get())
            size.append(scale.get())
            top.destroy()
            window_status.append('close')
            graphic()
        def color_selection():
            a=askcolor()
            colors.delete(0,END)
            colors.insert(INSERT,a[1])
        def shape_selection():
            if shape_types.get()=='oval':
                shapes.delete(0,END)
                shapes.insert(INSERT,'rectangle')
        if window_status[len(window_status)-1]=='close':
            def close_window():
                top.destroy()
                window_status.append('close')
            top=Toplevel(win)
            top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
            top.title('Graph settings')
            Frame(top,width=550,height=300,bg='#f5f5f5').pack()
            top.resizable(0,0)
            shape_types=StringVar()
            color_types=StringVar()
            size_values=StringVar()
            def font_selection():
                def accept_font(event):
                    font_types.delete(0,END)
                    index=font_list.get(font_list.curselection())
                    font_types.insert(INSERT,index)
                    graph_font.append(index)
                    font_list.destroy()
                font_list=Listbox(top,width=15,height=4,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                font_list.place(x=350,y=70)
                font_list.insert(END,'Arial','Calibri','Century Gothic','Courier New','Garamond', 'Helvetica','Lucida Sans','Tahoma','Times New Roman','Trebuchet MS','Verdana')
                font_list.bind("<Double-Button-1>",accept_font)
            def font_size_selection():
                def accept_size(event):
                    font_size.delete(0,END)
                    index=size_list.get(size_list.curselection())
                    font_size.insert(INSERT,index)
                    graph_font_size.append(index)
                    size_list.destroy()
                size_list=Listbox(top,width=15,height=4,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                size_list.place(x=350,y=110)
                size_list.insert(END,'6','8','9','10','11','12','14','16')
                size_list.bind("<Double-Button-1>",accept_size)
            Frame(top,width=250,height=200,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=10,y=15)
            Label(top,text='Main settings',font=font_type,bg='#f5f5f5',relief='ridge').place(x=15,y=5)
            Label(top,text='Color:',font=font_type,bg='#fcfcfc',relief='flat').place(x=20,y=35)
            Label(top,text='Shape:',font=font_type,bg='#fcfcfc',relief='flat').place(x=20,y=75)
            Label(top,text='Size:',font=font_type,bg='#fcfcfc',relief='flat').place(x=20,y=115)
            Frame(top,width=250,height=200,bg='#fcfcfc',borderwidth=2,relief='groove').place(x=280,y=15)
            Label(top,text='Labels',font=font_type,bg='#f5f5f5',relief='ridge').place(x=285,y=5)
            Radiobutton(top,text='None',font=font_type,bg='#fcfcfc',variable=label_display,value=0,relief='flat',cursor='hand2').place(x=290,y=35)
            Radiobutton(top,text='Visible',font=font_type,bg='#fcfcfc',variable=label_display,value=1,relief='flat',cursor='hand2').place(x=400,y=35)
            Label(top,text='Font:',font=font_type,bg='#fcfcfc',relief='flat').place(x=290,y=75)
            Label(top,text='Size:',font=font_type,bg='#fcfcfc',relief='flat').place(x=290,y=115)
            font_types=Entry(top,width=15,font=font_type,bg='#ffffff',relief='groove')
            font_types.place(x=350,y=70)
            font_types.insert(INSERT,'Trebuchet MS')
            font_size=Entry(top,width=15,font=font_type,bg='#ffffff',relief='groove')
            font_size.place(x=350,y=110)
            font_size.insert(INSERT,'6')
            font_button=Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=font_selection)
            font_button.place(x=450,y=65)
            size_button=Button(top,image=foto25,bg='#fcfcfc',relief='flat',overrelief='groove',cursor='hand2',command=font_size_selection)
            size_button.place(x=450,y=105)
            colors=Entry(top,width=15,font=font_type,bg='#ffffff',textvariable=color_types,relief='groove')
            colors.place(x=70,y=35)
            colors.insert(INSERT,'#8080c0')
            shapes=Entry(top,width=15,font=font_type,bg='#ffffff',textvariable=shape_types,relief='groove')
            shapes.place(x=70,y=75)
            shapes.insert(INSERT,'oval')
            scale=Scale(top,width=10,length=154,from_=0,to=30,tick=5,font=font_type,bg='#fcfcfc',cursor='hand2',orient='horizontal')
            scale.place(x=70,y=110)
            scale.set(15)
            Button(top,text='Select',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=color_selection).place(x=180,y=30)
            Button(top,text='Switch',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=shape_selection).place(x=180,y=70)
            Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=manuals).place(x=10,y=260)
            Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=accept_options).place(x=410,y=260)
            Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=cancel_graphic_options).place(x=470,y=260)
            window_status.append('open')
            top.protocol("WM_DELETE_WINDOW",close_window)
    def open_GeneQuery_project():
        try:
            file=askopenfilename(filetypes = [('GeneQuery', '.as')])
            aff.append(file)
            lefi=file[::-1]
            n1=lefi.find('/')
            lefi_name=lefi[0:n1]
            real_file_name=lefi_name[::-1]
            information_file=open('C:/PG/processes/data/projects/'+real_file_name,'r')
            content_information_file=information_file.readlines()
            genomic_technology.append(content_information_file[0][0:len(content_information_file[0])-1])
            active_workflow.append(content_information_file[0][0:len(content_information_file[0])-1])
            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                del(active_workflow[0])
                active_workflow.append('Affymetyrix')
            for i in range(1,len(content_information_file),1):
                analysis_type.append(content_information_file[i][0:len(content_information_file[i])-1])
            for i in range(1,len(content_information_file),1):
                workflow.append(content_information_file[i][0:len(content_information_file[i])-1])
            try:
                file_name.append(file)
                file_type.append('txt')
                message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                message_label.place(x=10,y=550)
                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                information_text.config(state=NORMAL)
                if genomic_technology[len(genomic_technology)-1]=='Custom':
                    information_text.insert(INSERT,'- Custom data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                    information_text.insert(INSERT,'- Custom data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Illumina':
                    information_text.insert(INSERT,'- Illumina data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                    information_text.insert(INSERT,'- Agilent data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Genepix':
                    information_text.insert(INSERT,'- Genepix data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Imagene':
                    information_text.insert(INSERT,'- Imagene data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Scanarray':
                    information_text.insert(INSERT,'- Scanarray data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Quantarray':
                    information_text.insert(INSERT,'- Quantarray data project imported.\n')
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    information_text.insert(INSERT,'- Affymetrix data project imported.\n')
                information_text.config(state=DISABLED)
                win.update()
                r.source('c:\\PG\\srs.R')
                warning=askquestion(title='Open Data',message='Are your data in Log Scale?')
                if warning=='yes':
                    r.projectcustom(file,'TRUE','txt',genomic_technology[len(genomic_technology)-1])
                    log_scale.append('TRUE')
                if warning=='no':
                    r.projectcustom(file,'FALSE','txt',genomic_technology[len(genomic_technology)-1])
                    log_scale.append('FALSE')
                content=open(file,'r')
                infiles=open('c:\\PG\\processes\\data\\result.txt','w')
                contents=content.readlines()
                column_number=(contents[0].count('\t'))+1
                infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
                infile=r.matrix(infile)
                for i in range(0,len(contents),1):
                    infiles.write(contents[i])
                infiles.close()
                for i in range(0,len(infile),1):
                    columns.append(infile[i][0][0])
                for i in range(1,len(infile[1][0]),1):
                    if log_scale[len(log_scale)-1]=='FALSE':
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                    if log_scale[len(log_scale)-1]=='TRUE':
                        if genomic_technology[len(genomic_technology)-1]=='Custom':
                            gene_id.append(str(infile[0][0][i]));data1.append(e**float(infile[1][0][i]));data2.append(e**float(infile[2][0][i]))
                        if genomic_technology[len(genomic_technology)-1]!='Custom':
                            if genomic_technology[len(genomic_technology)-1]!='qPCR' and genomic_technology[len(genomic_technology)-1]!='Agilent_CGH':
                                gene_id.append(str(infile[0][0][i]));data1.append(2**float(infile[1][0][i]));data2.append(2**float(infile[2][0][i]))
                            if genomic_technology[len(genomic_technology)-1]=='qPCR' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                                gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))                            
                win.update()
                for i in range(0,len(data1),1):
                    if data1[i]>values[len(values)-1]:
                        values.append(data1[i])
                    if data2[i]>values2[len(values2)-1]:
                        values2.append(data2[i])
                    if values2[len(values2)-1]>values[len(values)-1]:
                        values.append(values2[len(values2)-1])
                for i in range(0,len(data1),1):
                    new_data1.append((data1[i]*650)/values[len(values)-1])
                    new_data2.append((data2[i]*400)/values[len(values)-1])
                analysis_type.append('raw')
                data_table()
                win.update()
                sample1.insert(INSERT,columns[2])
                sample2.insert(INSERT,columns[1])
                message_label.destroy()
                graphic()
                workflow.append('preprocessing')
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    analysis_type.append('Affymetrix')
                    a=technology_class()
                    a.affymetrix()
                if genomic_technology[len(genomic_technology)-1]=='Illumina':
                    analysis_type.append('Illumina')
                    a=technology_class()
                    a.illumina()
                if genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Imagene' or genomic_technology[len(genomic_technology)-1]=='SMD':
                    analysis_type.append('Two channel')
                    a=technology_class()
                    a.two_channel()
                if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                    analysis_type.append('CGH')
                    tidy_up_mark.append('TRUE')
                    a=technology_class()
                    a.two_channel()
                if genomic_technology[len(genomic_technology)-1]=='Custom':
                    a=technology_class()
                    a.custom()
                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                    analysis_type.append('qPCR')
                    a=technology_class()
                    a.qPCR()
                if genomic_technology[len(genomic_technology)-1]!='Custom':
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                    if analysis_type.count('imputation')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if analysis_type.count('normalization')>=1:
                        if genomic_technology[len(genomic_technology)-1]!='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)   
                    if analysis_type.count('limma')>=1:
                        if genomic_technology[len(genomic_technology)-1]!='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                    if analysis_type.count('annotations')>=1:
                        if genomic_technology[len(genomic_technology)-1]!='qPCR':
                            Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                    if analysis_type.count('cluster')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=255)
                    if analysis_type.count('tidy_up')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    if analysis_type.count('segmentation')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                if genomic_technology[len(genomic_technology)-1]=='Custom':
                    Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=105)
                    if analysis_type.count('normalization')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=135)
                    if analysis_type.count('limma')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=165)
                    if analysis_type.count('annotations')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=195)
                    if analysis_type.count('cluster')>=1:
                        Label(win,image=foto17,bg='#fcfcfc',relief='flat').place(x=215,y=225)
                project.append('open')
            except:
                showwarning(title='Open file',message='Data format not valid. Importing process aborted.')
                new_project()
            r.source('c:\\PG\\srs.R')
            if genomic_technology[len(genomic_technology)-1]!='Affymetrix' and genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                r.project(file,'')
            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                r.project(file,real_file_name)
        except:
            pass
    def export_project():
        enrichment_results=IntVar()
        def cancel_project_export():
            top.destroy()
            window_status.append('close')
        def accept_project_export():
            top.destroy()
            top.update()
            window_status.append('close')
            if export_option.get()==0:
                file=asksaveasfilename(filetypes = [('GeneQuery', '.as')])
                infile=open(file,'w')
                files=open('c:\\PG\\processes\\data\\result.txt','r')
                lefi=file[::-1]
                n1=lefi.find('/')
                name=lefi[0:n1]
                file_name=name[::-1]
                information_file=open('C:/PG/processes/data/projects/'+file_name,'w')
                information_file.write(genomic_technology[len(genomic_technology)-1]+'\n')
                if analysis_type.count('normalization')>=1:
                    information_file.write('normalization'+'\n')
                if analysis_type.count('tidy_up')>=1:
                    information_file.write('tidy_up'+'\n')
                if analysis_type.count('segmentation')>=1:
                    information_file.write('segmentation'+'\n')
                if analysis_type.count('limma')>=1:
                    information_file.write('limma'+'\n')
                if analysis_type.count('sam')>=1:
                    information_file.write('sam'+'\n')
                if analysis_type.count('ttest')>=1:
                    information_file.write('ttest'+'\n')
                if analysis_type.count('cluster')>=1:
                    information_file.write('cluster'+'\n')
                if analysis_type.count('annotations')>=1:
                    information_file.write('annotations'+'\n')
                content=files.readlines()
                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                    hybridization_file=open('c:/PG/processes/data/projects/'+file_name+'_hyb.txt','w')
                    hyb_file=open('c:/PG/processes/data/hyb_result.txt','r')
                    hyb_content=hyb_file.readlines()
                    for i in range(0,len(hyb_content),1):
                        hybridization_file.write(hyb_content[i])
                    hybridization_file.close()     
                    hyb_file.close()
                for i in range (0,len(content),1):
                    infile.write(content[i])
                infile.close()
            if export_option.get()==1:
                try:
                    os.system('start excel c:\\PG\\processes\\data\\result.txt')
                    if enrichment_results.get()==1:
                        os.system('start excel c:\\PG\\processes\\data\\enrichment_result.txt')
                except:
                    warning=showinfo(title='Export Project',message='Impossible to open connection with Microsoft Excel')
            if export_option.get()==2:
                save_file()
                if enrichment_results.get()==1:
                    warning=showinfo(title='Export Project',message='Common saving process does not allow to save Biological Enrichment results')
        if window_status[len(window_status)-1]=='close' and genomic_technology[len(genomic_technology)-1]!='none':
            def close_window():
                top.destroy()
                window_status.append('close')
            top=Toplevel(win)
            top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
            top.title('Export Project')
            top.resizable(0,0)
            top.maxsize(300,250)
            export_option=IntVar()
            Frame(top,bg='#f5f5f5',width=300,height=250,relief='groove',borderwidth=2).pack()
            Frame(top,bg='#fcfcfc',width=270,height=170,relief='groove',borderwidth=2).place(x=10,y=20)
            Label(top,text='Data table',font=font_type,bg='#fcfcfc',relief='ridge').place(x=15,y=10)
            Radiobutton(top,text='Save as GeneQuery project',font=font_type,bg='#fcfcfc',variable=export_option,value=0,cursor='hand2').place(x=30,y=40)
            Radiobutton(top,text='Export to MS Excel',font=font_type,bg='#fcfcfc',variable=export_option,value=1,cursor='hand2').place(x=30,y=70)
            Radiobutton(top,text='Save data table as (*.txt, *.csv)',font=font_type,bg='#fcfcfc',variable=export_option,value=2,cursor='hand2').place(x=30,y=100)
            Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=manuals).place(x=10,y=210)
            Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=accept_project_export).place(x=170,y=210)
            Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',cursor='hand2',relief='groove',command=cancel_project_export).place(x=230,y=210) 
            window_status.append('open')
            top.protocol("WM_DELETE_WINDOW", close_window)
        if genomic_technology[len(genomic_technology)-1]=='none':
            warning=showinfo(title='Export Project',message='Dataset must be imported before exporting')
    def clipboard():
        question=askquestion('Copy','This action will copy all table content to clipboard\nDo you want to continue?')
        if question=='yes':
            s=Tk()
            s.withdraw()
            s.clipboard_clear()
            infile=open('C:\\PG\\processes\\data\\result.txt','r')
            content=infile.readlines()
            for i in range (0,len(content),1):
                s.clipboard_append(content[i])
            s.destroy()
    def clipboard_callback(event):
        clipboard()
    def open_file():
        if len(file_name)>0:
            warning=showinfo(title='Open Data',message='Save and close this project before opening a new one')
        if len(file_name)==0:
            del(file_name[1:len(file_name)])
            del(columns[1:len(columns)])
            file=askopenfilename(filetypes = [('all files', '.*'), ('text files', '.txt'),('csv files', '.csv')])
            lefi=file[::-1]
            if lefi[0:3]=='txt' or lefi[0:2]=='sa':
                try:
                    file_name.append(file)
                    file_type.append('txt')
                    message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                    message_label.place(x=10,y=550)
                    process.create_line(0,7,50,7,fill='#8080c0',width=15)
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- Custom data imported. Normalization is recommended if it was not previously applied.\n')
                    information_text.config(state=DISABLED)
                    win.update()
                    r.source('c:\\PG\\srs.R')
                    warning=askquestion(title='Open Data',message='Are your data in Log Scale?')
                    if warning=='yes':
                        r.custom(file,'TRUE','txt')
                        log_scale.append('TRUE')
                    if warning=='no':
                        r.custom(file,'FALSE','txt')
                        log_scale.append('FALSE')
                    content=open(file,'r')
                    infiles=open('c:\\PG\\processes\\data\\result.txt','w')
                    contents=content.readlines()
                    column_number=(contents[0].count('\t'))+1
                    infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
                    infile=r.matrix(infile)
                    for i in range(0,len(contents),1):
                        infiles.write(contents[i])
                    infiles.close()
                    for i in range(0,len(infile),1):
                        columns.append(infile[i][0][0])
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                    win.update()
                    for i in range(0,len(data1),1):
                        if data1[i]>values[len(values)-1]:
                            values.append(data1[i])
                        if data2[i]>values2[len(values2)-1]:
                            values2.append(data2[i])
                        if values2[len(values2)-1]>values[len(values)-1]:
                            values.append(values2[len(values2)-1])
                    for i in range(0,len(data1),1):
                        new_data1.append((data1[i]*650)/values[len(values)-1])
                        new_data2.append((data2[i]*400)/values[len(values)-1])
                    if active_workflow[len(active_workflow)-1]!='none':
                        presentation_text.destroy()
                        a=technology_class()
                        a.custom()
                    genomic_technology.append('Custom')
                    analysis_type.append('raw')
                    data_table()
                    win.update()
                    sample1.insert(INSERT,columns[2])
                    sample2.insert(INSERT,columns[1])
                    message_label.destroy()
                    graphic()
                except:
                    showwarning(title='Open file',message='Data format not valid. Importing process aborted.')
                    new_project()
            if lefi[0:3]=='vsc':
                try:
                    file_name.append(file)
                    file_type.append('vsc')                                     
                    message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                    message_label.place(x=10,y=550)
                    process.create_line(0,7,50,7,fill='#8080c0',width=15)
                    information_text.config(state=NORMAL)
                    information_text.insert(INSERT,'- Custom data imported. Normalization is recommended if it was not previously applied.\n')
                    information_text.config(state=DISABLED)
                    win.update()
                    r.source('c:\\PG\\srs.R')
                    content=open(file,'r')
                    infiles=open('c:\\PG\\processes\\data\\result.txt','w')
                    contents=content.readlines()
                    column_number=(contents[0].count(','))+1
                    infile=r.read_table(file,sep=',',col_names=col_name[0:column_number])
                    infile=r.matrix(infile)
                    for i in range(0,len(contents),1):
                        index=contents[i].replace(',','\t')
                        infiles.write(index)
                    infiles.close()
                    for i in range(0,len(infile),1):
                        columns.append(infile[i][0][0])
                    for i in range(1,len(infile[1][0]),1):
                        gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                    win.update()
                    for i in range(0,len(data1),1):
                        if data1[i]>values[len(values)-1]:
                            values.append(data1[i])
                        if data2[i]>values2[len(values2)-1]:
                            values2.append(data2[i])
                        if values2[len(values2)-1]>values[len(values)-1]:
                            values.append(values2[len(values2)-1])
                    for i in range(0,len(data1),1):
                        new_data1.append((data1[i]*650)/values[len(values)-1])
                        new_data2.append((data2[i]*400)/values[len(values)-1])
                    if active_workflow[len(active_workflow)-1]!='none':
                        presentation_text.destroy()
                        a=technology_class()
                        a.custom()
                    warning=askquestion(title='Open Data',message='Are your data in Log Scale?')
                    if warning=='yes':
                        r.custom(file,'TRUE','vsc')
                        log_scale.append('TRUE')
                    if warning=='no':
                        r.custom(file,'FALSE','vsc')
                        log_scale.append('FALSE')
                    genomic_technology.append('Custom')
                    analysis_type.append('raw')
                    data_table()
                    win.update()
                    sample1.insert(INSERT,columns[2])
                    sample2.insert(INSERT,columns[1])
                    message_label.destroy()
                    graphic()
                except:
                    showwarning(title='Open file',message='Data format not valid. Importing process aborted.')
                    new_project()
    def paste_clipboard():
        question=askquestion('Paste','This action will copy all clipboard content to GeneQuery\nDo you want to continue?')
        s=Tk()
        contenido=s.clipboard_get()
        if len(file_name)>0:
            warning=showinfo(title='Data Pivoting',message='Save and close this project before pasting')
        if question=='yes' and len(file_name)==0:
            infile=open('C:\\PG\\processes\\data\\clipboard_result.txt','w')
            for i in range (0,len(contenido),1):
                infile.write(contenido[i])
            infile.close()
            s.destroy()
            try:
                del(file_name[1:len(file_name)])
                del(columns[0:len(columns)])
                file='C:\\PG\\processes\\data\\clipboard_result.txt'
                file_name.append(file)
                file_type.append('txt')
                message_label=Label(win,text='Importing Data...',font=font_type3,bg='#fcfcfc',relief='flat')
                message_label.place(x=10,y=550)
                process.create_line(0,7,50,7,fill='#8080c0',width=15)
                information_text.config(state=NORMAL)
                information_text.insert(INSERT,'- Custom data pasted from clipboard.\n')
                information_text.config(state=DISABLED)
                win.update()
                r.source('c:\\PG\\srs.R')
                warning=askquestion(title='Open Data',message='Are your data in Log Scale?')
                if warning=='yes':
                    r.custom(file,'TRUE','txt')
                    log_scale.append('TRUE')
                if warning=='no':
                    r.custom(file,'FALSE','txt')
                    log_scale.append('FALSE')
                content=open(file,'r')
                infiles=open('c:\\PG\\processes\\data\\result.txt','w')
                contents=content.readlines()
                column_number=(contents[0].count('\t'))+1
                infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
                infile=r.matrix(infile)
                for i in range(0,len(contents),1):
                    infiles.write(contents[i])
                infiles.close()
                for i in range(0,len(infile),1):
                    columns.append(infile[i][0][0])
                for i in range(1,len(infile[1][0]),1):
                    gene_id.append(str(infile[0][0][i]));data1.append(float(infile[1][0][i]));data2.append(float(infile[2][0][i]))
                win.update()
                for i in range(0,len(data1),1):
                    if data1[i]>values[len(values)-1]:
                        values.append(data1[i])
                    if data2[i]>values2[len(values2)-1]:
                        values2.append(data2[i])
                    if values2[len(values2)-1]>values[len(values)-1]:
                        values.append(values2[len(values2)-1])
                for i in range(0,len(data1),1):
                    new_data1.append((data1[i]*650)/values[len(values)-1])
                    new_data2.append((data2[i]*400)/values[len(values)-1])
                if active_workflow[len(active_workflow)-1]!='none':
                    presentation_text.destroy()
                    a=technology_class()
                    a.custom()
                genomic_technology.append('Custom')
                analysis_type.append('raw')
                data_table()
                win.update()
                sample1.insert(INSERT,columns[2])
                sample2.insert(INSERT,columns[1])
                message_label.destroy()
                graphic()
            except:
                warning=showwarning(title='Paste',message='Impossible to paste clipboard in GeneQuery')
                message_label.destroy()
                process.destroy()
                try:
                    os.remove('c:\\PG\\processes\\data\\clipboard_result.txt')
                except:
                    pass
        if question=='no':
            pass
    def paste_callback(event):
        paste_clipboard()
    def save_file():
        if genomic_technology[len(genomic_technology)-1]!='none':
            if analysis_type[len(analysis_type)-1]=='none' or analysis_type[len(analysis_type)-1]=='raw':
                try:
                    file=asksaveasfilename(filetypes = [('all files', '.*'), ('text files', '.txt'),('csv files', '.csv')])
                    files=open(file_name[len(file_name)-1],'r')
                    lefi=file[::-1]
                    content=files.readlines()
                    infile=open(file,'w')
                    if lefi[0:3]!='vsc':
                        for i in range(0,len(content),1): 
                            infile.write(content[i])
                        infile.close()
                    if lefi[0:3]=='vsc':
                        for i in range(0,len(content),1):
                            index=content[i].replace('\t',',')
                            infile.write(index)
                        infile.close()
                except:
                    warning=showwarning(title='Save Data',message='Data saving process aborted')
            if analysis_type[len(analysis_type)-1]!='none'and analysis_type[len(analysis_type)-1]!='raw':
                try:
                    file=asksaveasfilename(filetypes = [('all files', '.*'), ('text files', '.txt'),('csv files', '.csv')])
                    files=open('c:\\PG\\processes\\data\\result.txt','r')
                    lefi=file[::-1]
                    content=files.readlines()
                    infile=open(file,'w')
                    if lefi[0:3]!='vsc':
                        for i in range (0,len(content),1):
                            infile.write(content[i])
                        infile.close()
                    if lefi[0:3]=='vsc':
                        for i in range (0,len(content),1):
                            index=content[i].replace('\t',',')
                            infile.write(index)
                        infile.close()
                except:
                    warning=showwarning(title='Save Data',message='Data saving process aborted')
        if genomic_technology[len(genomic_technology)-1]=='none':
            warning=showinfo(title='Save Data',message='No data to be exported or saved')
    def genomic_analysis():
        microarrays_list=['close']
        if len(file_name)>0:
            warning=showinfo(title='Open Data',message='Save and close this project before opening a new one')
        if len(file_name)==0:
            if window_status[len(window_status)-1]=='close':
                def close_window():
                    top.destroy()
                    window_status.append('close')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.maxsize(500,500)
                top.title('Analysis Workflow')
                Frame(top,width=600,height=600,bg='#f5f5f5').pack()
                technology_chosen=StringVar()
                def cancel_analysis_workflow():
                    top.destroy()
                    window_status.append('close')
                def technology_selection():
                    if microarray_list[len(microarray_list)-1]!='qPCR' :
                        def technology_selected(event):
                            index=technology_list.get(technology_list.curselection())
                            technology_entry.config(state=NORMAL,disabledbackground='#ffffff')
                            technology_entry.delete(0,END)
                            technology_entry.insert(INSERT,index)
                            technology_entry.config(state=DISABLED,disabledbackground='#ffffff')
                            technology_list.destroy()
                            microarrays_list.append('close')
                            if genomic_technology[len(genomic_technology)-1]=='qPCR':
                                genomic_analysis()
                        if microarrays_list[len(microarrays_list)-1]=='close':
                            microarrays_list.append('open')
                            technology_list=Listbox(top,width=40,height=8,bg='#ffffff',font=font_type,cursor='hand2',relief='groove')
                            technology_list.place(x=125,y=10)
                            technology_list.insert(END,'Affymetrix','Agilent','Custom','Genepix','Illumina','Imagene','Quantarray','Scanarray')
                            technology_list.bind("<Double-Button-1>",technology_selected)
                def import_technology():
                    top.destroy()
                    win.update()
                    def accept_importing():
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                           index=int(missing_ct.get())
                           qpcr_missing.append(index)
                        if genomic_technology[len(genomic_technology)-1]=='Affymetrix':
                            if len(folder_entry.get())>=1:
                                infile=open('c:\\PG\\processes\\data\\Affymetrix\\design.txt','r')
                                content=infile.readlines()
                                if len(content)>=2:
                                    active_workflow.append('Affymetrix')
                                    up.destroy()
                                    up.update()
                                    a=technology_class()
                                    a.affymetrix()
                                    analysis_type.append('Affymetrix')
                                if len(content)<=1:
                                    warning=showinfo(title='Import files',message='Select two files at least before starting Affymetrix data analysis')
                        if genomic_technology[len(genomic_technology)-1]=='Illumina':
                            if len(folder_entry.get())>=1:
                                infile=open('c:\\PG\\processes\\data\\Illumina\\design.txt','r')
                                content=infile.readlines()
                                if len(content)>=1:
                                    active_workflow.append('Illumina')
                                    up.destroy()
                                    up.update()
                                    a=technology_class()
                                    a.illumina()
                                    analysis_type.append('Illumina')
                                if len(content)<1:
                                    warning=showinfo(title='Import files',message='Select one file at least before starting Illumina data analysis')
                        if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray'or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                            if len(folder_entry.get())>=1:
                                infile=open('c:\\PG\\processes\\data\\twochannel\\design.txt','r')
                                content=infile.readlines()
                            if len(content)>1:
                                    active_workflow.append('Two channel')
                                    up.destroy()
                                    up.update()
                                    a=technology_class()
                                    a.two_channel()
                                    analysis_type.append('Two channel')
                            if len(content)<=1:
                                    warning=showinfo(title='Import files',message='Select one file at least before starting two channel data analysis')
                        if genomic_technology[len(genomic_technology)-1]=='Agilent_CGH':
                            if len(folder_entry.get())>=1:
                                infile=open('c:\\PG\\processes\\data\\twochannel\\design.txt','r')
                                content=infile.readlines()
                            if len(content)>1:
                                    active_workflow.append('CGH')
                                    up.destroy()
                                    up.update()
                                    a=technology_class()
                                    a.two_channel()
                                    analysis_type.append('CGH')
                            if len(content)<=1:
                                    warning=showinfo(title='Import files',message='Select one file at least before starting two channel data analysis')
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            if len(folder_entry.get())>=1:
                                infile=open('c:\\PG\\processes\\data\\qPCR\\design.txt','r')
                                content=infile.readlines()
                                if len(content)>=1:
                                    active_workflow.append('qPCR')
                                    up.destroy()
                                    up.update()
                                    a=technology_class()
                                    a.qPCR()
                                    analysis_type.append('qPCR')
                                if len(content)<1:
                                    warning=showinfo(title='Import files',message='Select one file at least before starting qPCR data analysis')
                        try:
                            file_name.append(folder[len(folder)-1])
                        except:
                            warning=showinfo(title='Import files',message='Select one file at least before starting the analysis')
                    def cancel_analysis_workflow2():
                        up.destroy()
                        del(genomic_technology[0:len(genomic_technology)-1])
                        genomic_technology.append('none')
                        up.update()
                    def browse_folder():
                        folder_entry.config(state=NORMAL,disabledbackground='#ffffff')
                        folder_entry.delete(0,END)
                        if browse_status[len(browse_status)-1]=='close':
                            browse_status.append('open')
                            if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                file=askopenfilenames(filetypes = [('cel', '.cel'),('chp', '.chp')])
                            if genomic_technology[len(genomic_technology)-1]!='Affymetrix' and genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                                if genomic_technology[len(genomic_technology)-1]=='Genepix':
                                    file=askopenfilenames(filetypes = [('Genepix files','.gpr'),('text files', '.txt'),('csv','.csv')])
                                if genomic_technology[len(genomic_technology)-1]!='Genepix' and genomic_technology[len(genomic_technology)-1]!='Imagene':
                                    file=askopenfilenames(filetypes = [('text files', '.txt'),('csv','.csv')])
                                if genomic_technology[len(genomic_technology)-1]=='Imagene':
                                    file=askopenfilenames(filetypes = [('Imagene files', '.cy*')])
                            if len(file)==0:
                                browse_status.append('close')
                            try:
                                if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                    files=open('c:\\PG\\processes\\data\\Affymetrix\\design.txt','w')
                                if genomic_technology[len(genomic_technology)-1]=='Illumina':
                                    files=open('c:\\PG\\processes\\data\\Illumina\\design.txt','w')
                                if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                                    files=open('c:\\PG\\processes\\data\\Twochannel\\design.txt','w')
                                    files.write('Samples\tFileName\n')
                                    files.close()
                                    files=open('c:\\PG\\processes\\data\\Twochannel\\design.txt','a')
                                if genomic_technology[len(genomic_technology)-1]=='qPCR':
                                    files=open('c:\\PG\\processes\\data\\qPCR\\design.txt','w')
                                if genomic_technology[len(genomic_technology)-1]!='Custom':
                                    infile=file[0][::-1]
                                    n1=infile.find('/')
                                    content=infile[n1:len(infile)]
                                    folder_entry.insert(INSERT,content[::-1])
                                    folder_entry.config(state=DISABLED,disabledbackground='#ffffff')
                                    folder.append(content[::-1])
                                    for i in range(0,len(file),1):
                                        file1.insert(END,file[i])
                                        infiles=file[i][::-1]
                                        n1=infiles.find('/')
                                        contents=infiles[4:n1]
                                        files.write(contents[::-1]+'\t'+file[i]+'\n')
                                    files.close()
                            except:
                                top.destroy()
                    def clear_list():
                        file1.delete(0,END)
                        if genomic_technology[len(genomic_technology)-1]=='Affymetrix' or genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                            files=open('c:\\PG\\processes\\data\\Affymetrix\\design.txt','w')
                        if genomic_technology[len(genomic_technology)-1]=='Illumina':
                            files=open('c:\\PG\\processes\\data\\Illumina\\design.txt','w')
                        if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='AgilentCGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                                files=open('c:\\PG\\processes\\data\\Twochannel\\design.txt','w')
                                files.write('Samples\tFileName\n')
                                files.close()
                        if genomic_technology[len(genomic_technology)-1]=='qPCR':
                            files=open('c:\\PG\\processes\\data\\qPCR\\design.txt','w')
                    up=Toplevel(win)
                    up.resizable (0,0)
                    up.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                    def folder_button_label(event):
                        def clear_folder_label(event):
                            folder_label.destroy()
                        folder_label=Label(up,text='Browse dataset',font=font_type,bg='#fcfcfc',relief='groove')
                        folder_label.place(x=380,y=35)
                        folder_button.bind("<Leave>",clear_folder_label)
                    up.title('Analysis '+ genomic_technology[len(genomic_technology)-1])
                    Frame(up,width=500,height=500,bg='#f5f5f5').pack()
                    Label(up,text='Selected folder:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=10)
                    Label(up,text='Files to analyze:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=65)
                    folder_entry=Entry(up,width=40,bg='#ffffff',font=font_type,relief='groove')
                    folder_entry.place(x=120,y=10)
                    folder_entry.config(state=DISABLED,disabledbackground='#ffffff')
                    folder_button=Button(up,image=foto1,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=browse_folder)
                    folder_button.place(x=365,y=5)
                    folder_button.bind("<Enter>",folder_button_label)
                    Button(up,text='Help',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=manuals).place(x=25,y=450)
                    Button(up,text='Clear',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=clear_list).place(x=90,y=450)
                    Button(up,text='Continue',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=accept_importing).place(x=360,y=450)
                    Button(up,text='Cancel',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=cancel_analysis_workflow2).place(x=425,y=450)
                    file1=Listbox(up,width=75,height=15,bg='#ffffff',font=font_type,relief='groove',cursor='hand2')
                    file1.place(x=25,y=90)
                    if genomic_technology[len(genomic_technology)-1]=='Affymetrix':                    
                        Label(up,text='Summarization method:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=380)                    
                        Radiobutton(up,text='MAS 5.0',font=font_type,bg='#f5f5f5',variable=analysis_method,value=0,relief='flat',cursor='hand2').place(x=170,y=380)                    
                        Radiobutton(up,text='RMA',font=font_type,bg='#f5f5f5',variable=analysis_method,value=1,relief='flat',cursor='hand2').place(x=270,y=380)
                    if genomic_technology[len(genomic_technology)-1]=='Illumina':
                        def browse_controls():
                            try:
                                file=askopenfilename()
                                control_probes.delete(0,END)
                                control_probes.insert(INSERT,file)
                                illumina_controls_file.append(control_probes.get())
                            except:
                                control_probes.delete(0,END)
                                control_probes.insert(INSERT,'none')
                        def control_folder_label(event):
                            def clear_control_label(event):
                                control_label.destroy()
                            control_label=Label(up,text='Browse control file',font=font_type,bg='#fcfcfc',relief='groove')
                            control_label.place(x=390,y=405)
                            control_button.bind("<Leave>",clear_control_label)
                        Label(up,text='Control Probe file:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=380)
                        control_probes=Entry(up,width=40,font=font_type,bg='#fcfcfc',relief='groove')
                        control_probes.place(x=130,y=380)
                        control_probes.config(state=NORMAL,disabledbackground='#fcfcfc')
                        control_probes.insert(INSERT,'none')
                        control_probes.config(state=DISABLED,disabledbackground='#fcfcfc')
                        control_button=Button(up,image=foto1,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                        control_button.place(x=375,y=375)
                        control_button.bind("<Enter>",control_folder_label)
                    if genomic_technology[len(genomic_technology)-1]=='Genepix' or genomic_technology[len(genomic_technology)-1]=='Quantarray' or genomic_technology[len(genomic_technology)-1]=='Scanarray' or genomic_technology[len(genomic_technology)-1]=='Agilent' or genomic_technology[len(genomic_technology)-1]=='Agilent_CGH' or genomic_technology[len(genomic_technology)-1]=='Imagene':
                        def background_methods(event):
                            def background_selection(event):
                                index=background_list.get(background_list.curselection())
                                background_list.destroy()
                                background_method.config(state=NORMAL,disabledbackground='#ffffff')
                                background_method.delete(0,END)
                                background_method.insert(INSERT,index)
                                background_method.config(state=DISABLED,disabledbackground='#ffffff')
                                backgrounds_list.append ('close')
                            if backgrounds_list[len(backgrounds_list)-1]=='close':
                                backgrounds_list.append('open')
                                background_list=Listbox(up,height=4,width=40,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                                background_list.place(x=150,y=360)
                                background_list.insert(END,'none','subtract','minimum','normexp')
                                background_list.bind("<Double-Button-1>",background_selection)
                        def normalization_methods(event):
                            def within_normalization_selection(event):
                                index=normalization_list.get(normalization_list.curselection())
                                normalization_list.destroy()
                                normalization_method.config(state=NORMAL,disabledbackground='#ffffff')
                                normalization_method.delete(0,END)
                                normalization_method.insert(INSERT,index)
                                normalization_method.config(state=DISABLED,disabledbackground='#ffffff')
                                normalizations_list.append('close')
                            if normalizations_list[len(normalizations_list)-1]=='close':
                                normalizations_list.append('open')
                                normalization_list=Listbox(up,height=3,width=40,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                                normalization_list.place(x=150,y=390)
                                normalization_list.insert(END,'none','median','loess')
                                normalization_list.bind("<Double-Button-1>",within_normalization_selection)
                        Label(up,text='Background correction:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=360)
                        background_method=Entry(up,width=40,font=font_type,bg='#fcfcfc',textvariable=background,relief='groove')
                        background_method.place(x=150,y=360)
                        background_method.config(state=NORMAL,disabledbackground='#ffffff')
                        background_method.delete(0,END)
                        background_method.insert(INSERT,'none')
                        background_method.config(state=DISABLED,disabledbackground='#ffffff')
                        background_button=Button(up,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                        background_button.place(x=395,y=355)
                        background_button.bind("<Button-1>",background_methods)
                        Label(up,text='Within normalization:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=390)
                        normalization_method=Entry(up,width=40,font=font_type,bg='#fcfcfc',textvariable=within_normalization,relief='groove')
                        normalization_method.place(x=150,y=390)
                        normalization_method.config(state=NORMAL,disabledbackground='#ffffff')
                        normalization_method.delete(0,END)
                        normalization_method.insert(INSERT,'none')
                        normalization_method.config(state=DISABLED,disabledbackground='#ffffff')
                        normalization_button=Button(up,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                        normalization_button.place(x=395,y=385)
                        normalization_button.bind("<Button-1>",normalization_methods)
                    if genomic_technology[len(genomic_technology)-1]=='qPCR':
                        def aggregate_replicates():
                            def aggregate_selection(event):
                                index=aggregation_methods.get(aggregation_methods.curselection())
                                aggregation_methods.destroy()
                                aggregation.config(state=NORMAL,disabledbackground='#ffffff')
                                aggregation.delete(0,END)
                                aggregation.insert(INSERT,index)
                                aggregation.config(state=DISABLED,disabledbackground='#ffffff')
                            aggregation_methods=Listbox(up,width=10,height=2,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                            aggregation_methods.place(x=180,y=390)
                            aggregation_methods.insert(END,'median','mean')
                            aggregation_methods.bind("<Double-Button-1>", aggregate_selection)
                        Label(up,text='Change missing Ct values by:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=360)
                        missing_values=Entry(up,width=10,font=font_type,bg='#fcfcfc',textvariable=missing_ct,relief='groove',cursor='hand2')
                        missing_values.place(x=180,y=360)
                        missing_values.insert(INSERT,'35')
                        Label(up,text='Replicates aggregation:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=390)
                        aggregation=Entry(up,width=10,font=font_type,bg='#fcfcfc',relief='groove')
                        aggregation.place(x=180,y=390)
                        aggregation.config(state=NORMAL,disabledbackground='#ffffff')
                        aggregation.insert(INSERT,'median')
                        aggregation.config(state=DISABLED,disabledbackground='#ffffff')
                        aggregate_method.append(aggregation.get())
                        Button(up,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=aggregate_replicates).place(x=245,y=385)
                def accept_technology():
                    genomic_technology.append(technology_chosen.get())
                    window_status.append('close')
                    if genomic_technology[len(genomic_technology)-1]=='Custom':
                        active_workflow.append('Custom')
                        top.destroy()
                        open_file()
                    if genomic_technology[len(genomic_technology)-1]!='Custom':
                        import_technology()
                Label(top,text='Select technology:',font=font_type,bg='#f5f5f5',relief='flat').place(x=25,y=10)
                technology_entry=Entry(top,width=40,bg='#ffffff',font=font_type,textvariable=technology_chosen,relief='groove')
                technology_entry.place(x=125,y=10)
                technology_entry.config(state=NORMAL,disabledbackground='#ffffff')
                technology_entry.insert(INSERT,microarray_list[len(microarray_list)-1])
                technology_entry.config(state=DISABLED,disabledbackground='#ffffff')
                if microarray_list[len(microarray_list)-1]!='Agilent_CGH':
                    technology_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=technology_selection)
                if microarray_list[len(microarray_list)-1]=='Agilent_CGH':
                    technology_button=Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2')
                technology_button.place(x=370,y=5)
                if microarray_list[len(microarray_list)-1]!='Citometry':
                    Label(top,text='Gene Expression Data Analysis guides the user throughout a step-by-step workflow.',font=font_type,bg='#f5f5f5',fg='#8080c0',relief='flat').place(x=25,y=50)
                    Label(top,text='Click on each step to get more detailed information.',font=font_type,bg='#f5f5f5',fg='#8080c0',relief='flat').place(x=25,y=70)
                Frame(top,width=450,height=300,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=25,y=100)
                if microarray_list[len(microarray_list)-1]!='qPCR' and microarray_list[len(microarray_list)-1]!='Citometry':
                    Label(top,text='1. Microarrays Preprocessing.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=120)
                    Label(top,text='Uploads and preprocess raw data.', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=250,y=125)
                    Label(top,text='2. Quality Control.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=170)
                    Label(top,text='Metrics displayed to help users monitoring the quality of', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=170,y=175)
                    Label(top,text='raw data, examine the replicates and the signal distribution of all arrays.',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=195)
                    Label(top,text='3. Data Normalization.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=220)
                    Label(top,text='Removes systematic errors from data, making ', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=195,y=225)
                    Label(top,text='Microarrays comparable among them. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=245)
                    if microarray_list[len(microarray_list)-1]!='Agilent_CGH':
                        Label(top,text='4. Differential Expression Analysis.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=270)
                        Label(top,text='Compares Microarray datas to find', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=285,y=275)
                        Label(top,text='genes differentially expressed between two different conditions. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=295)
                    if microarray_list[len(microarray_list)-1]=='Agilent_CGH':
                        Label(top,text='4. Tide Up.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=270)
                        Label(top,text='Reorders microarrays data information associating chromosomes', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=125,y=275)
                        Label(top,text=' and possitions to each clone in the analysis. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=295)
                    if microarray_list[len(microarray_list)-1]!='Agilent_CGH':                       
                        Label(top,text='5. Functional Analysis.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=320)
                        Label(top,text='Provides information of all biological processes', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=195,y=325)
                        Label(top,text='associated to those genes differentially expressed. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=345)
                    if microarray_list[len(microarray_list)-1]=='Agilent_CGH':                       
                        Label(top,text='5. Segmentation.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=320)
                        Label(top,text='Summarizes the replicates, reorders the data and divides', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=165,y=325)
                        Label(top,text='the genome in the corresponding regions to estimate alterations. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=345)                                            
                if microarray_list[len(microarray_list)-1]=='qPCR':
                    Label(top,text='1. Import Ct values.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=120)
                    Label(top,text='Uploads and preprocess raw Ct values', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=180,y=125)
                    Label(top,text='2. Quality Control.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=170)
                    Label(top,text='Metrics displayed to help users monitoring the quality of', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=170,y=175)
                    Label(top,text='Ct values, examine the replicates and the signal distribution of all plates.',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=195)
                    Label(top,text='3. Data Imputation.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=220)
                    Label(top,text='Replaces the missing Ct values for a gene that forms', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=175,y=225)
                    Label(top,text='part of the same sample group . ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=245)
                    Label(top,text='4. Data Normalization.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=270)
                    Label(top,text='Removes systematic errors from data, making', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=200,y=275)
                    Label(top,text='plates comparable among them. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=295)                    
                    Label(top,text='5. Differential Expression Analysis.',bg='#fcfcfc',fg='#8080c0',font=font_type3,relief='flat').place(x=35,y=320)
                    Label(top,text='Compares different plates to find', bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=290,y=325)
                    Label(top,text='genes differentially expressed between two different conditions. ',bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat').place(x=55,y=345)
                Button(top,text='Help',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=manuals).place(x=25,y=450)
                Button(top,text='Continue',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=accept_technology).place(x=360,y=450)
                Button(top,text='Cancel',bg='#fcfcfc',width=7,font=font_type,relief='groove',cursor='hand2',command=cancel_analysis_workflow).place(x=425,y=450)
                window_status.append('open')
                top.protocol("WM_DELETE_WINDOW",close_window)
    def exit_genomics_office():
        warning=askquestion(title='Exit GeneQuery',message='Do you want to close and leave the current project?')
        if warning=='yes':
            try:
                os.remove('c:\\PG\\processes\\data\\raw_data.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\boxplot.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\boxplot_max.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\cor_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\normalized_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\normalized_a_values.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\hyb_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\control_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\affymetrix_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\affymetrix_normalized_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\enrichment_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\illumina_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\custom_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\box_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\twochannel_normalized_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\annotation.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\report.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\main_graph.ps')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\quality_graph.ps')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\clipboard_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\pivot_result.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\design.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\samples.txt')
            except:
                pass
            try:
                os.remove('c:\\PG\\processes\\data\\project_data.txt')
            except:
                pass
            win.destroy()
    def cluster_data():
        if len(file_name)>0:
            a=preprocess_class()
            a.hierarchical_clustering()
    def normalize_data():
        if len(file_name)>0:
            a=preprocess_class()
            a.data_normalization()
    def impute_data():
        if len(file_name)>0:
            a=preprocess_class()
            a.data_imputation()
    def summarize_data():
        if len(file_name)>0:
            a=preprocess_class()
            a.data_summarization()
    def transpose_data():
        if len(file_name)>0:
            a=preprocess_class()
            a.data_transposition()
    def limma_data_analysis():
        if len(file_name)>0:
            a=statistic_class()
            a.limma()
    def sam_data_analysis():
        if len(file_name)>0:
            a=statistic_class()
            a.sam()
    def ttest_data_analysis():
        if len(file_name)>0:
            a=statistic_class()
            a.t_test()
    def data_enrichment():
        if len(file_name)>0:
            a=functional_analysis()
            a.biological_enrichment()
    def annotation_retriever():
        if len(file_name)>0:
            a=functional_analysis()
            a.retrieve_annotations()
    def load_custom_script():
            a=preprocess_class()
            a.custom_script()
    def open_file_label(event):
        def clear_open_label(event):
            open_label.destroy()
        open_label=Label(text='Open',font=font_type,bg='#fcfcfc',relief='groove')
        open_label.place(x=15,y=35)
        open_files.bind("<Leave>",clear_open_label)
    def save_file_label(event):
        def clear_save_label(event):
            save_label.destroy()
        save_label=Label(text='Save',font=font_type,bg='#fcfcfc',relief='groove')
        save_label.place(x=45,y=35)
        save_files.bind("<Leave>",clear_save_label)
    def new_file_label(event):
        def clear_new_label(event):
            new_label.destroy()
        new_label=Label(text='Clear',font=font_type,bg='#fcfcfc',relief='groove')
        new_label.place(x=80,y=35)
        new.bind("<Leave>",clear_new_label)
    def printer_file_label(event):
        def clear_printer_label(event):
            printer_label.destroy()
        printer_label=Label(text='Copy',font=font_type,bg='#fcfcfc',relief='groove')
        printer_label.place(x=120,y=35)
        printer.bind("<Leave>",clear_printer_label)
    def export_file_label(event):
        def clear_export_label(event):
            export_label.destroy()
        export_label=Label(text='Export',font=font_type,bg='#fcfcfc',relief='groove')
        export_label.place(x=150,y=35)
        export.bind("<Leave>",clear_export_label)
    def genomics_file_label(event):
        def clear_genomics_label(event):
            genomics_label.destroy()
        genomics_label=Label(text='Microarrays',font=font_type,bg='#fcfcfc',relief='groove')
        genomics_label.place(x=190,y=35)
        genomics.bind("<Leave>",clear_genomics_label)
    def tools_file_label(event):
        def clear_tools_label(event):
            tools_label.destroy()
        tools_label=Label(text='Tools',font=font_type,bg='#fcfcfc',relief='groove')
        tools_label.place(x=250,y=35)
        tools.bind("<Leave>",clear_tools_label)
    def report_file_label(event):
        def clear_report_label(event):
            report_label.destroy()
        report_label=Label(text='Custom Scripts',font=font_type,bg='#fcfcfc',relief='groove')
        report_label.place(x=310,y=35)
        report.bind("<Leave>",clear_report_label)
    def table_file_label(event):
        def clear_table_label(event):
            table_label.destroy()
        table_label=Label(text='User Workflow',font=font_type,bg='#fcfcfc',relief='groove')
        table_label.place(x=220,y=35)
        table.bind("<Leave>",clear_table_label)
    def create_report_label(event):
        def clear_report_label(event):
            report_label.destroy()
        report_label=Label(text='Create report',font=font_type,bg='#fcfcfc',relief='groove')
        report_label.place(x=280,y=35)
        create_report.bind("<Leave>",clear_report_label)
    def question_file_label(event):
        def clear_question_label(event):
            question_label.destroy()
        question_label=Label(text='Help',font=font_type,bg='#fcfcfc',relief='groove')
        question_label.place(x=345,y=35)
        question.bind("<Leave>",clear_question_label)
    def manuals():
        try:
            os.system('start AcroRd32 C:\PG\Text\Help\help_document.pdf')
        except:
            pass
    def GeneQuery_tools():
        def quit_tools():
            top.destroy()
            window_status.append('close')
        if window_status[len(window_status)-1]=='close' and genomic_technology[len(genomic_technology)-1]!='none':
            def clustering(event):
                quit_tools()
                cluster_data()
            def normalize(event):
                quit_tools()
                normalize_data()
            def impute(event):
                quit_tools()
                impute_data()
            def summarize(event):
                quit_tools()
                summarize_data()
            def transpose(event):
                quit_tools()
                transpose_data()
            def limma_analysis(event):
                quit_tools()
                limma_data_analysis()
            def sam_analysis(event):
                quit_tools()
                sam_data_analysis()
            def ttest_analysis(event):
                quit_tools()
                ttest_data_analysis()
            def enrichment(event):
                quit_tools()
                data_enrichment()
            def ID_annotation(event):
                quit_tools()
                annotation_retriever()
            window_status.append('open')
            top=Toplevel(win)
            top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
            top.resizable(0,0)
            top.title('GeneQuery Tools')
            Frame(top,width=330,height=440,bg='#f5f5f5',relief='groove',borderwidth=2).pack()
            Frame(top,width=315,height=110,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=15)
            Frame(top,width=315,height=110,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=135)
            Frame(top,width=315,height=110,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=255)
            Label(top,image=foto13,bg='#fcfcfc',cursor='hand2').place(x=10,y=20)
            Label(top,image=foto33,bg='#fcfcfc',cursor='hand2').place(x=10,y=140)
            Label(top,image=foto34,bg='#fcfcfc',cursor='hand2').place(x=10,y=260)
            Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=manuals).place(x=10,y=390)
            Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=quit_tools).place(x=265,y=390)
            Label(top,text='Data Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=20)
            normalization=Label(top,text='- Normalization',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            normalization.place(x=65,y=95)
            normalization.bind("<Button-1>",normalize)
            imputation=Label(top,text='- Imputation',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            imputation.place(x=65,y=70)
            imputation.bind("<Button-1>",impute)
            summarization=Label(top,text='- Summarization',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            summarization.place(x=65,y=45)
            summarization.bind("<Button-1>",summarize)
            transposition=Label(top,text='- Transposition',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            transposition.place(x=170,y=45)
            transposition.bind("<Button-1>",transpose)
            Label(top,text='Statistical analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=140)
            limma=Label(top,text='- Linear Models for Microarrays (LIMMA)',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            limma.place(x=65,y=165)
            limma.bind("<Button-1>",limma_analysis)
            sam=Label(top,text='- Significant Analysis of Microarrays (SAM)',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            sam.place(x=65,y=190)
            sam.bind("<Button-1>",sam_analysis)
            ttest=Label(top,text='- Unpaired T-Test',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            ttest.place(x=65,y=215)
            ttest.bind("<Button-1>",ttest_analysis)
            Label(top,text='Cluster and Functional analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=260)
            cluster=Label(top,text='- Hierarchical Clustering',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            cluster.place(x=65,y=285)
            cluster.bind("<Button-1>",clustering)
            annotation=Label(top,text='- Retrieve annotations',font=font_type1,bg='#fcfcfc',fg='#8080c0',cursor='hand2')
            annotation.place(x=65,y=310)
            annotation.bind("<Button-1>",ID_annotation)
            top.protocol("WM_DELETE_WINDOW",quit_tools)
        if genomic_technology[len(genomic_technology)-1]=='none':
            warning=showinfo(title='GeneQuery Tools',message='Dataset must be imported before displaying GeneQuery Tools')
    def user_workflow():
        position1=[75,95,115];position2=[170,190,210];position3=[255,275]
        def quit_user_workflow():
            top.destroy()
            window_status.append('close')
            user_workflow_status.append('none')
        def accept_user_workflow():
            if LIMMA.get()==1 and SAM.get()==1 and TTEST.get()==1:
                warning=showinfo(title='GeneQuery user workflow',message='Please select just one differential expression analysis method.')
            if LIMMA.get()==0 and SAM.get()==1 and TTEST.get()==1:
                warning=showinfo(title='GeneQuery user workflow',message='Please select just one differential expression analysis method.')
            if LIMMA.get()==1 and SAM.get()==0 and TTEST.get()==1:
                warning=showinfo(title='GeneQuery user workflow',message='Please select just one differential expression analysis method.')
            if LIMMA.get()==1 and SAM.get()==1 and TTEST.get()==0:
                warning=showinfo(title='GeneQuery user workflow',message='Please select just one differential expression analysis method.')
            if (LIMMA.get()==1 and SAM.get()==0 and TTEST.get()==0) or (LIMMA.get()==0 and SAM.get()==1 and TTEST.get()==0) or (LIMMA.get()==0 and SAM.get()==0 and TTEST.get()==1):
                presentation_text.destroy()
                win.update()
                def save_custom_workflow(event):
                    file=asksaveasfilename(filetypes = [('GeneQuery', '.as'),('all files', '.*')])
                    infile=open(file,'w')                   
                    infile.write(str(summarize.get())+'\n'+str(impute.get())+'\n'+str(normalize.get())+'\n'+str(LIMMA.get())+'\n'+str(SAM.get())+'\n'+str(TTEST.get())+'\n'+str(clustering.get())+'\n'+str(biological_enrichment.get()))
                    infile.close()
                def summarize_callback(event):
                    summarize_data()
                def impute_callback(event):
                    impute_data()
                def normalize_callback(event):
                    normalize_data()
                def limma_data_analysis_callback(event):
                    limma_data_analysis()
                def sam_data_analysis_callback(event):
                    sam_data_analysis()
                def ttest_data_analysis_callback(event):
                    ttest_data_analysis()
                def cluster_callback(event):
                    cluster_data()
                def enrichment_callback(event):
                    annotation_retriever()
                Frame1=Frame(win,width=247,height=100,bg='#f5f5f5',relief='groove',borderwidth=2)
                Frame1.place(x=7,y=40)
                Frame2=Frame(win,width=247,height=80,bg='#f5f5f5',relief='groove',borderwidth=2)
                Frame2.place(x=7,y=142)
                Frame3=Frame(win,width=247,height=80,bg='#f5f5f5',relief='groove',borderwidth=2)
                Frame3.place(x=7,y=225)
                Label(win,text='1. Data Preprocessing', font=font_type3,bg='#f5f5f5',fg='#8080c0').place(x=9,y=50)
                label1=Label(win,text='- Summarization', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label1.bind("<Button-1>",summarize_callback)
                label2=Label(win,text='- Imputation', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label2.bind("<Button-1>",impute_callback)
                label3=Label(win,text='- Normalization', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label3.bind("<Button-1>",normalize_callback)
                Label(win,text='2. Statistical Analysis', font=font_type3,bg='#f5f5f5',fg='#8080c0').place(x=9,y=145)
                label4=Label(win,text='- Linear Models for Microarrays (LIMMA)', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label4.bind("<Button-1>",limma_data_analysis_callback)
                label5=Label(win,text='- Significant Analysis of Microarrays (SAM)', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label5.bind("<Button-1>",sam_data_analysis_callback)
                label6=Label(win,text='- Unpaired T-Test', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label6.bind("<Button-1>",ttest_data_analysis_callback)
                Label(win,text='3. Cluster and Functional analysis', font=font_type3,bg='#f5f5f5',fg='#8080c0').place(x=9,y=230)
                label7=Label(win,text='- Hierarchical Clustering', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label7.bind("<Button-1>",cluster_callback)
                label8=Label(win,text='- Annotation retriever', font=font_type,bg='#f5f5f5',fg='#8080c0',cursor='hand2')
                label8.bind("<Button-1>",enrichment_callback)
                save_workflow=Label(win,text='Save template?',font=font_type,bg='#ffffff',fg='#8080c0',cursor='hand2')
                save_workflow.place(x=165, y=305)
                save_workflow.bind("<Button-1>",save_custom_workflow)
                if summarize.get()==1:
                    label1.place(x=15,y=position1[0])
                    del(position1[0])
                if impute.get()==1:
                    label2.place(x=15,y=position1[0])
                    del(position1[0])
                if normalize.get()==1:
                    label3.place(x=15,y=position1[0])
                    del(position1[0])
                if LIMMA.get()==1 and SAM.get()==0 and TTEST.get()==0:
                    label4.place(x=15,y=position2[0])
                    del(position2[0])
                if SAM.get()==1 and LIMMA.get()==0 and TTEST.get()==0:
                    label5.place(x=15,y=position2[0])
                    del(position2[0])
                if TTEST.get()==1 and SAM.get()==0 and LIMMA.get()==0:
                    label6.place(x=15,y=position2[0])
                    del(position2[0])
                if clustering.get()==1:
                    label7.place(x=15,y=position3[0])
                    del(position3[0])
                if biological_enrichment.get()==1:
                    label8.place(x=15,y=position3[0])
                    del(position3[0])
                top.destroy()
                window_status.append('close')
        if window_status[len(window_status)-1]=='close' and genomic_technology[len(genomic_technology)-1]!='none' and user_workflow_status[len(user_workflow_status)-1]=='none':
            if active_workflow[len(active_workflow)-1]=='none':
                user_workflow_status.append('open')
                top=Toplevel(win)
                top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
                top.resizable(0,0)
                top.title('GeneQuery User Workflow')
                window_status.append('open')
                def open_template():
                    try:
                        file=askopenfilename(filetypes = [('GeneQuery', '.as'),('all files', '.*')])
                        infile=open(file,'r')
                        content=infile.readlines()
                        summarize.set(int(content[0]))
                        impute.set(int(content[1]))
                        normalize.set(int(content[2]))
                        LIMMA.set(int(content[3]))
                        SAM.set(int(content[4]))
                        TTEST.set(int(content[5]))
                        clustering.set(int(content[6]))
                        biological_enrichment.set(int(content[7]))
                        user_workflow()
                    except:
                        pass
                Frame(top,width=340,height=460,bg='#f5f5f5',relief='groove',borderwidth=2).pack()
                Frame(top,width=325,height=120,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=15)
                Frame(top,width=325,height=120,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=145)
                Frame(top,width=325,height=120,bg='#fcfcfc',relief='groove',borderwidth=2).place(x=5,y=275)
                Label(top,image=foto13,bg='#fcfcfc',cursor='hand2').place(x=10,y=20)
                Label(top,image=foto33,bg='#fcfcfc',cursor='hand2').place(x=10,y=150)
                Label(top,image=foto34,bg='#fcfcfc',cursor='hand2').place(x=10,y=280)
                Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=manuals).place(x=10,y=410)
                Button(top,text='Template',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=open_template).place(x=75,y=410)
                Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=accept_user_workflow).place(x=210,y=410)
                Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=quit_user_workflow).place(x=275,y=410)
                Label(top,text='Data Preprocessing',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=20)
                summarization=Checkbutton(top,text='- Summarization',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=summarize,relief='flat',cursor='hand2')
                summarization.place(x=65,y=45)
                imputation=Checkbutton(top,text='- Imputation',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=impute,relief='flat',cursor='hand2')
                imputation.place(x=65,y=70)
                normalization=Checkbutton(top,text='- Normalization',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=normalize,relief='flat',cursor='hand2')
                normalization.place(x=65,y=95)
                Label(top,text='Statistical analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=150)
                limma=Checkbutton(top,text='- Linear Models for Microarrays (LIMMA)',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=LIMMA,relief='flat',cursor='hand2')
                limma.place(x=65,y=175)
                sam=Checkbutton(top,text='- Significant Analysis of Microarrays (SAM)',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=SAM,relief='flat',cursor='hand2')
                sam.place(x=65,y=200)
                ttest=Checkbutton(top,text='- Unpaired T-Test',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=TTEST,relief='flat',cursor='hand2')
                ttest.place(x=65,y=225)
                Label(top,text='Cluster and Functional analysis',font=font_type3,bg='#fcfcfc',fg='#8080c0',cursor='hand2').place(x=45,y=280)
                cluster=Checkbutton(top,text='- Hierarchical Clustering',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=clustering,relief='flat',cursor='hand2')
                cluster.place(x=65,y=305)
                enrichment=Checkbutton(top,text='- Annotation retriever',font=font_type1,bg='#fcfcfc',fg='#8080c0',variable=biological_enrichment,relief='flat',cursor='hand2')
                enrichment.place(x=65,y=330)
                top.protocol("WM_DELETE_WINDOW",quit_user_workflow)
    #Disposicion botones en la barra de herramientas y Logotipo
    open_files=Button(image=foto1,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=open_file)
    open_files.place(x=5,y=5)
    open_files.bind("<Enter>",open_file_label)
    save_files=Button(image=foto2,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=save_file)
    save_files.place(x=35,y=5)
    save_files.bind("<Enter>",save_file_label)
    new=Button(image=foto15,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=start_new_project)
    new.place(x=65,y=5)
    new.bind("<Enter>",new_file_label)
    Label(text='|',bg='#c9c9f1',relief='flat',font=font_type).place(x=95,y=8)
    printer=Button(image=foto5,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=clipboard)
    printer.place(x=105,y=5)
    printer.bind("<Enter>",printer_file_label)
    export=Button(image=foto4,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=export_project)
    export.place(x=135,y=5)
    export.bind("<Enter>",export_file_label)
    Label(text='|',bg='#c9c9f1',relief='flat',font=font_type).place(x=165,y=8)
    genomics=Button(image=foto8,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=genomic_analysis)
    genomics.place(x=175,y=5)
    genomics.bind("<Enter>",genomics_file_label)
    tools=Button(image=foto9,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=GeneQuery_tools)
    tools.place(x=235,y=5)
    tools.bind("<Enter>",tools_file_label)
    create_report=Button(image=foto35,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=reporting)
    create_report.place(x=265,y=5)
    create_report.bind("<Enter>", create_report_label)
    report=Button(image=foto10,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=load_custom_script)
    report.place(x=295,y=5)
    report.bind("<Enter>",report_file_label)
    table=Button(image=foto11,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=user_workflow)
    table.place(x=205,y=5)
    table.bind("<Enter>",table_file_label)
    Label(text='|',bg='#c9c9f1',relief='flat',font=font_type).place(x=325,y=8)
    question=Button(image=foto16,bg='#c9c9f1',activebackground='#c9c9f1',relief='flat',overrelief='groove',cursor='hand2',command=manuals)
    question.place(x=335,y=5)
    question.bind("<Enter>",question_file_label)
    options_panel=Frame(win,width=250,height=660,bg='#fcfcfc',relief='groove',borderwidth=2)
    options_panel.place(x=5,y=35)
    if width_resolution[len(width_resolution)-1]>1024:
        result_panel=Frame(win,width=750+x_resolution,height=660,bg='#ffffff',relief='groove',borderwidth=2)
    if width_resolution[len(width_resolution)-1]<=1024:
        result_panel=Frame(win,width=750,height=660,bg='#ffffff',relief='groove',borderwidth=2)
    result_panel.place(x=265,y=35)
    process=Canvas(win,width=200,height=14,bg='#ffffff',relief='flat',bd=-2)
    process.place(x=10,y=590)
    Label(win,text='----------------- Information panel -----------------',bg='#fcfcfc',font=font_type,relief='flat').place(x=10,y=320)
    information_text=Text(win,width=35,height=11,bg='#fcfcfc',fg='#8080c0',font=font_type,relief='flat', wrap='word')
    information_text.place(x=10,y=340)
    information_text.config(state=DISABLED)
    scrollbar=Scrollbar(win,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
    scrollbar.place(x=225,y=340,height=170)
    information_text.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=information_text.yview)
    presentation_text=Text(win,width=39,height=29,font=font_type1,bg='#fcfcfc',fg='#8080c0',wrap='word',relief='flat')
    presentation_text.place(x=10,y=45)
    presentation_text.image_create(END, image=foto28, align=CENTER)
    text=open('c:\\PG\\text\\presentation.txt','r')
    content=text.readlines()
    for i in range (0,len(content),1):
        presentation_text.insert(INSERT,content[i])
    Label(text='Intensity of',font=font_type,bg='#ffffff',relief='flat').place(x=320+x_resolution/2,y=220)
    presentation_text.config(state=DISABLED)
    def sample1_display(event):
        if analysis_type[len(analysis_type)-1]!='none':
            def display1_selection(event):
                sample1.config(state=NORMAL,disabledbackground='#ffffff')
                sample1.delete(0,END)
                index=display1.get(display1.curselection())
                if index!='GeneID' and index!='gene_id':
                    sample1.insert(INSERT,index)
                    display1.destroy()
                    graph_status=Canvas(win,width=100,height=10,bg='#fcfcfc',borderwidth=1,relief='flat')
                    graph_status.place(x=580+x_resolution/2,y=430)
                    graph_status.create_line(0,6,25,6,fill='#8080c0',width=10)
                    win.update()
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                    infile=r.matrix(infile)
                    del(data1[0:len(data1)]);del(data2[0:len(data2)])
                    del(new_data1[0:len(new_data1)]);del(new_data2[0:len(new_data2)])
                    del(values[0:len(values)]);del(values2[0:len(values2)-1])
                    values.append(0);values2.append(0)
                    graph_status.create_line(0,6,50,6,fill='#8080c0',width=10)
                    win.update()
                    for i in range(0,len(infile),1):
                        if sample1.get()==(infile[i][0][0]):
                            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    data1.append(2**float(infile[i][0][j]))
                            if genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    if str(infile[i][0][j])=='NA':
                                        infile[i][0][j]=0
                                    if str(infile[i][0][j])!='NA':
                                        data1.append(float(infile[i][0][j]))
                        if sample2.get()==(infile[i][0][0]):
                            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    data2.append(2**float(infile[i][0][j]))
                            if genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    if str(infile[i][0][j])=='NA':
                                        infile[i][0][j]=0
                                    if str(infile[i][0][j])!='NA':
                                        data2.append(float(infile[i][0][j]))
                    for i in range(0,len(data1),1):
                        if data1[i]>values[len(values)-1]:
                            values.append(data1[i])
                        if data2[i]>values2[len(values2)-1]:
                            values2.append(data2[i])
                    if values2[len(values2)-1]>values[len(values)-1]:
                        values.append(values2[len(values2)-1])
                    graph_status.create_line(0,6,75,6,fill='#8080c0',width=10)
                    win.update()
                    for i in range(0,len(data1),1):
                        new_data1.append((data2[i]*650)/values[len(values)-1])
                        new_data2.append((data1[i]*400)/values[len(values)-1])
                    graph_status.create_line(0,6,100,6,fill='#8080c0',width=10)
                    win.update()
                    graph_status.destroy()
                    graphic()
                    sample1.config(state=DISABLED,disabledbackground='#fcfcfc')
                if index=='GeneID' or index=='gene_id':
                    warning=showinfo(title='Scatter Plot',message='Elements represented within the scatter plot should be numeric')
            display1=Listbox(win,width=20,height=6,font=font_type,bg='#ffffff',relief='groove',cursor='hand2')
            display1.place(x=390+x_resolution/2,y=220)
            for i in range(0,len(columns),1):
                display1.insert(END,columns[i])
            display1.bind("<Double-Button-1>",display1_selection)
    def sample2_display(event):
        if analysis_type[len(analysis_type)-1]!='none':
            def display2_selection(event):
                sample2.config(state=NORMAL,disabledbackground='#ffffff')
                sample2.delete(0,END)
                index=display2.get(display2.curselection())
                if index!='GeneID' and index!='gene_id':
                    sample2.config(state=NORMAL,disabledbackground='#ffffff')
                    sample2.insert(INSERT,index)
                    display2.destroy()
                    graph_status=Canvas(win,width=100,height=10,bg='#fcfcfc',borderwidth=1,relief='flat')
                    graph_status.place(x=580+x_resolution/2,y=430)
                    graph_status.create_line(0,6,25,6,fill='#8080c0',width=10)
                    win.update()
                    infile=r.read_table('c:\\PG\\processes\\data\\result.txt',sep='\t')
                    infile=r.matrix(infile)
                    graph_status.create_line(0,6,50,6,fill='#8080c0',width=10)
                    del(data1[0:len(data1)]);del(data2[0:len(data2)])
                    del(new_data1[0:len(new_data1)]);del(new_data2[0:len(new_data2)])
                    del(values[0:len(values)]);del(values2[0:len(values2)-1])
                    values.append(0);values2.append(0)
                    for i in range(0,len(infile),1):
                        if sample1.get()==(infile[i][0][0]):
                            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':    
                                for j in range(1,len(infile[1][0]),1):
                                    data1.append(2**float(infile[i][0][j]))
                            if genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    if str(infile[i][0][j])=='NA':
                                        infile[i][0][j]=0
                                    if str(infile[i][0][j])!='NA':
                                        data1.append(float(infile[i][0][j]))
                        if sample2.get()==(infile[i][0][0]):
                            if genomic_technology[len(genomic_technology)-1]=='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    data2.append(2**float(infile[i][0][j]))
                            if genomic_technology[len(genomic_technology)-1]!='All-Affymetrix':
                                for j in range(1,len(infile[1][0]),1):
                                    if str(infile[i][0][j])=='NA':
                                        infile[i][0][j]=0
                                    if str(infile[i][0][j])!='NA':
                                        data2.append(float(infile[i][0][j]))
                    for i in range(0,len(data1),1):
                        if data1[i]>values[len(values)-1]:
                            values.append(data1[i])
                        if data2[i]>values2[len(values2)-1]:
                            values2.append(data2[i])
                    if values2[len(values2)-1]>values[len(values)-1]:
                        values.append(values2[len(values2)-1])
                    graph_status.create_line(0,6,75,6,fill='#8080c0',width=10)
                    win.update()
                    for i in range(0,len(data1),1):
                        new_data1.append((data2[i]*650)/values[len(values)-1])
                        new_data2.append((data1[i]*400)/values[len(values)-1])
                    graph_status.create_line(0,6,100,6,fill='#8080c0',width=10)
                    win.update()
                    graph_status.destroy()
                    graphic()
                    sample2.config(state=DISABLED,disabledbackground='#fcfcfc')
                if index=='GeneID' or index=='gene_id':
                    warning=showinfo(title='Scatter Plot',message='Elements represented within the scatter plot should be numeric')
            display2=Listbox(win,width=20,height=6,font=font_type,bg='#ffffff',relief='groove',cursor='hand2')
            display2.place(x=590+x_resolution/2,y=575)
            for i in range(0,len(columns),1):
                display2.insert(END,columns[i])
            display2.bind("<Double-Button-1>",display2_selection)
    def graphic_position(event):
        position_x.append(event.x)
        position_y.append(event.y)
    sample1=Entry(win,width=20,bg='#ffffff',fg='#0000ff',font=font_type,relief='flat',cursor='hand2')
    sample1.place(x=390+x_resolution/2,y=220)
    sample1.bind("<Button-1>",sample1_display)
    Label(text='Intensity of',font=font_type,bg='#ffffff',relief='flat').place(x=500+x_resolution/2,y=665)
    sample2=Entry(win,width=20,bg='#ffffff',fg='#0000ff',font=font_type,relief='flat',cursor='hand2')
    sample2.place(x=590+x_resolution/2,y=665)
    sample2.bind("<Button-1>",sample2_display)   
    graph_panel=Canvas(win,width=650,height=400,bg='#ffffff',cursor='crosshair',relief='groove',borderwidth=2)
    graph_panel.place(x=320+x_resolution/2,y=240)
    graph_panel.bind("<Button-3>",graphic_properties)
    graph_panel.bind("<Motion>",graphic_position)
    graph_panel.create_line(0,40,650,40,fill='#c5c5c5');graph_panel.create_line(0,80,650,80,fill='#c5c5c5')
    graph_panel.create_line(0,120,650,120,fill='#c5c5c5');graph_panel.create_line(0,160,650,160,fill='#c5c5c5')
    graph_panel.create_line(0,200,650,200,fill='#c5c5c5');graph_panel.create_line(0,240,650,240,fill='#c5c5c5')
    graph_panel.create_line(0,280,650,280,fill='#c5c5c5');graph_panel.create_line(0,320,650,320,fill='#c5c5c5')
    graph_panel.create_line(0,360,650,360,fill='#c5c5c5')
    graph_panel.create_line(65,0,65,400,fill='#c5c5c5');graph_panel.create_line(130,0,130,400,fill='#c5c5c5')
    graph_panel.create_line(195,0,195,400,fill='#c5c5c5');graph_panel.create_line(260,0,260,400,fill='#c5c5c5')
    graph_panel.create_line(325,0,325,400,fill='#c5c5c5');graph_panel.create_line(390,0,390,400,fill='#c5c5c5')
    graph_panel.create_line(455,0,455,400,fill='#c5c5c5');graph_panel.create_line(520,0,520,400,fill='#c5c5c5')
    graph_panel.create_line(585,0,585,400,fill='#c5c5c5')
    y0=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y0.place(x=280+x_resolution/2,y=230)
    y1=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y1.place(x=280+x_resolution/2,y=270)
    y2=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y2.place(x=280+x_resolution/2,y=310)
    y3=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y3.place(x=280+x_resolution/2,y=350)
    y4=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y4.place(x=280+x_resolution/2,y=390)
    y5=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y5.place(x=280+x_resolution/2,y=430)
    y6=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y6.place(x=280+x_resolution/2,y=470)
    y7=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y7.place(x=280+x_resolution/2,y=510)
    y8=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y8.place(x=280+x_resolution/2,y=550)
    y9=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y9.place(x=280+x_resolution/2,y=590)
    y10=Entry(win,width=5,bg='#ffffff',font=font_type,relief='flat')
    y10.place(x=280+x_resolution/2,y=630)
    x1=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x1.place(x=370+x_resolution/2,y=650)
    x2=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x2.place(x=435+x_resolution/2,y=650)
    x3=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x3.place(x=500+x_resolution/2,y=650)
    x4=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x4.place(x=565+x_resolution/2,y=650)
    x5=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x5.place(x=630+x_resolution/2,y=650)
    x6=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x6.place(x=695+x_resolution/2,y=650)
    x7=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x7.place(x=760+x_resolution/2,y=650)
    x8=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x8.place(x=825+x_resolution/2,y=650)
    x9=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x9.place(x=890+x_resolution/2,y=650)
    x10=Entry(win,width=6,bg='#ffffff',font=font_type,relief='flat')
    x10.place(x=955+x_resolution/2,y=650)
    x1.config(state=DISABLED,disabledbackground='#ffffff');x2.config(state=DISABLED,disabledbackground='#ffffff');x3.config(state=DISABLED,disabledbackground='#ffffff')
    x4.config(state=DISABLED,disabledbackground='#ffffff');x5.config(state=DISABLED,disabledbackground='#ffffff');x6.config(state=DISABLED,disabledbackground='#ffffff')
    x7.config(state=DISABLED,disabledbackground='#ffffff');x8.config(state=DISABLED,disabledbackground='#ffffff');x9.config(state=DISABLED,disabledbackground='#ffffff');x10.config(state=DISABLED,disabledbackground='#ffffff')
    y1.config(state=DISABLED,disabledbackground='#ffffff');y2.config(state=DISABLED,disabledbackground='#ffffff');y3.config(state=DISABLED,disabledbackground='#ffffff')
    y4.config(state=DISABLED,disabledbackground='#ffffff');y5.config(state=DISABLED,disabledbackground='#ffffff');y6.config(state=DISABLED,disabledbackground='#ffffff')
    y7.config(state=DISABLED,disabledbackground='#ffffff');y8.config(state=DISABLED,disabledbackground='#ffffff');y9.config(state=DISABLED,disabledbackground='#ffffff');y0.config(state=DISABLED,disabledbackground='#ffffff')
    data_table()
    def affymetrix_technology():
        microarray_list.append('Affymetrix')
        genomic_analysis()
    def illumina_technology():
        microarray_list.append('Illumina')
        genomic_analysis()
    def custom_technology():
        microarray_list.append('Custom')
        genomic_analysis()
    def genepix_technology():
        microarray_list.append('Genepix')
        genomic_analysis()
    def agilent_technology():
        microarray_list.append('Agilent')
        genomic_analysis()
    def CGH_technology():
        microarray_list.append('Agilent_CGH')
        genomic_analysis()
    def imagene_technology():
        microarray_list.append('Imagene')
        genomic_analysis()
    def scanarray_technology():
        microarray_list.append('Scanarray')
        genomic_analysis()
    def quantarray_technology():
        microarray_list.append('Quantarray')
        genomic_analysis()
    def qPCR_technology():
        microarray_list.append('qPCR')
        genomic_analysis()
    def user_case():
        try:
            os.system('start AcroRd32 C:\PG\Text\Help\user_case.pdf')
        except:
            pass
    def about_GeneQuery():
        if window_status[len(window_status)-1]=='close':
            def quit_about_GeneQuery():
                window_status.append('close')
                top.destroy()
            def contact_info(event):
                os.system('start iexplore http://about.me/joaquin.panadero')
            top=Toplevel(win)
            top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
            top.title('About GeneQuery')
            top.resizable(0,0)
            Frame(top,width=400,height=400,bg='#f5f5f5',borderwidth=2).pack()
            Label(top,image=foto28,bg='#fcfcfc',relief='ridge').place(x=15,y=20)
            Label(top,text='Version: 1.0.0',font=font_type,bg='#f5f5f5').place(x=300,y=70)
            Label(top,text='Product: GeneQuery for genomics', font=font_type,bg='#f5f5f5').place(x=10,y=210)
            Label(top,text='Edition: Basic', font=font_type,bg='#f5f5f5').place(x=10,y=230)
            Label(top,text='Packages: Affymetrix, Illumina, Two channel, qPCR, Annotation', font=font_type,bg='#f5f5f5').place(x=10,y=250)
            Label(top,text='Build date: 31/01/2011', font=font_type,bg='#f5f5f5').place(x=10,y=270)
            Label(top,text='For support please contact with:',font=font_type,bg='#f5f5f5').place(x=10,y=320)
            contact=Label(top,text='joaquin.panadero@hotmail.com',font=font_type,bg='#f5f5f5',fg='#0000ff',cursor='hand2')
            contact.place(x=180,y=320)
            contact.bind("<Button-1>",contact_info)
            Button(top,text='Close',font=font_type,width=7,bg='#fcfcfc',relief='groove',cursor='hand2',command=quit_about_GeneQuery).place(x=330,y=360)
            window_status.append('open')
            top.protocol("WM_DELETE_WINDOW",quit_about_GeneQuery)
    def pivot_data():
        def start_pivot():
            if len(columns)>0:
                genomic_technology.append('Custom')
                analysis_type.append('raw')
                top.destroy()
                a=preprocess_class()
                a.data_pivoting()
        def cancel_pivot():
            top.destroy()
        def open_long_format():
            file=askopenfilename()
            pivot_file_name.config(state=NORMAL,disabledbackground='#fcfcfc')
            pivot_file_name.insert(INSERT,file)
            pivot_file_name.config(state=DISABLED,disabledbackground='#fcfcfc')
            file_name.append(file)
            file_type.append('txt')
            win.update()
            r.source('c:\\PG\\srs.R')
            content=open(file,'r')
            infiles=open('c:\\PG\\processes\\data\\pivot_result.txt','w')
            contents=content.readlines()
            column_number=(contents[0].count('\t'))+1
            infile=r.read_table(file,sep='\t',col_names=col_name[0:column_number])
            infile=r.matrix(infile)
            for i in range(0,len(contents),1):
                infiles.write(contents[i])
            infiles.close()
            for i in range(0,len(infile),1):
                columns.append(infile[i][0][0])
            columns_name.insert(INSERT,columns[0])
            rows_name.insert(INSERT,columns[0])
            data_wide.insert(INSERT,columns[0])
        def aggregate_method():
            def aggregate_selection(event):
                index=aggregation_list.get(aggregation_list.curselection())
                aggregation_list.destroy()
                aggregation.config(state=NORMAL,disabledbackground='#fcfcfc')
                aggregation.delete(0,END)
                aggregation.insert(INSERT,index)
                aggregation.config(state=DISABLED,disabledbackground='#fcfcfc')
                aggregates_list.append('close')
            if aggregates_list[len(aggregates_list)-1]=='close':
                aggregation_list=Listbox(top,width=15,height=2,font=font_type,bg='#fcfcfc',relief='groove')
                aggregation_list.place(x=115,y=250)
                aggregation_list.insert(END,'median','mean')
                aggregation_list.bind("<Double-Button-1>", aggregate_selection)
                aggregates_list.append('open')
        def select_columns():
            if len(file_name)>=1:
                def column_selected(event):
                    index=columns_list.get(columns_list.curselection())
                    columns_list.destroy()
                    columns_scrollbar.destroy()
                    columns_name.config(state=NORMAL,disabledbackground='#fcfcfc')
                    columns_name.delete(0,END)
                    columns_name.insert(INSERT,index)
                    columns_name.config(state=DISABLED,disabledbackground='#fcfcfc')
                    pivot_columns.append(columns_name.get())
                    column_list.append('close')
                if column_list[len(column_list)-1]=='close':
                    columns_list=Listbox(top,width=40,height=5,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                    columns_list.place(x=45,y=95)
                    columns_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    columns_scrollbar.place(x=290,y=95,height=90)
                    columns_list.config(yscrollcommand=columns_scrollbar.set)
                    columns_scrollbar.config(command=columns_list.yview)
                    for i in range(0,len(columns),1):
                        columns_list.insert(END,columns[i])
                    columns_list.bind("<Double-Button-1>",column_selected)
                    column_list.append('open')
        def select_rows():
            if len(file_name)>=1:
                def row_selected(event):
                    index=rows_list.get(rows_list.curselection())
                    rows_list.destroy()
                    rows_scrollbar.destroy()
                    rows_name.config(state=NORMAL,disabledbackground='#fcfcfc')
                    rows_name.delete(0,END)
                    rows_name.insert(INSERT,index)
                    rows_name.config(state=DISABLED,disabledbackground='#fcfcfc')
                    pivot_rows.append(rows_name.get())
                    row_list.append('close')
                if row_list[len(row_list)-1]=='close':
                    rows_list=Listbox(top,width=40,height=5,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                    rows_list.place(x=45,y=155)
                    rows_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    rows_scrollbar.place(x=290,y=155,height=90)
                    rows_list.config(yscrollcommand=rows_scrollbar.set)
                    rows_scrollbar.config(command=rows_list.yview)
                    for i in range(0,len(columns),1):
                        rows_list.insert(END,columns[i])
                    rows_list.bind("<Double-Button-1>",row_selected)
                    row_list.append('open')
        def select_datas():
            if len(file_name)>=1:
                def data_selected(event):
                    index=datas_list.get(datas_list.curselection())
                    datas_list.destroy()
                    datas_scrollbar.destroy()
                    data_wide.config(state=NORMAL,disabledbackground='#fcfcfc')
                    data_wide.delete(0,END)
                    data_wide.insert(INSERT,index)
                    data_wide.config(state=DISABLED,disabledbackground='#fcfcfc')
                    pivot_datas.append(data_wide.get())
                    value_list.append('close')
                if value_list[len(value_list)-1]=='close':
                    datas_list=Listbox(top,width=40,height=5,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2')
                    datas_list.place(x=45,y=215)
                    datas_scrollbar=Scrollbar(top,relief='flat',activerelief='flat',troughcolor='#8080c0',elementborderwidth=1,cursor='hand2')
                    datas_scrollbar.place(x=290,y=215,height=90)
                    datas_list.config(yscrollcommand=datas_scrollbar.set)
                    datas_scrollbar.config(command=datas_list.yview)
                    for i in range(0,len(columns),1):
                        datas_list.insert(END,columns[i])
                    datas_list.bind("<Double-Button-1>",data_selected)
                    value_list.append('open')
        if len(file_name)==0:
            top=Toplevel(win)
            top.wm_iconbitmap('C:\\PG\\Icons\\GeneQuery.ico')
            top.title('Pivot data')
            top.resizable(0,0)
            Frame(top,width=350,height=370,bg='#f5f5f5',borderwidth=2).pack()
            Label(top,text='File:',font=font_type,bg='#f5f5f5').place(x=10,y=35)
            Label(top,image=foto38,bg='#f5f5f5',relief='flat').place(x=10,y=90)
            Label(top,image=foto36,bg='#f5f5f5',relief='flat').place(x=10,y=150)
            Label(top,image=foto37,bg='#f5f5f5',relief='flat').place(x=10,y=210)
            pivot_file_name=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
            pivot_file_name.place(x=45,y=35)
            pivot_file_name.config(state=DISABLED,disabledbackground='#fcfcfc')
            columns_name=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
            columns_name.place(x=45,y=95)
            columns_name.config(state=DISABLED,disabledbackground='#fcfcfc')
            rows_name=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
            rows_name.place(x=45,y=155)
            rows_name.config(state=DISABLED,disabledbackground='#fcfcfc')
            data_wide=Entry(top,width=40,font=font_type,bg='#fcfcfc',relief='groove')
            data_wide.place(x=45,y=215)
            data_wide.config(state=DISABLED,disabledbackground='#fcfcfc')
            Button(top,image=foto1,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=open_long_format).place(x=290,y=30)
            Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=select_columns).place(x=290,y=90)
            Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=select_rows).place(x=290,y=150)
            Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=select_datas).place(x=290,y=210)
            Label(top,text='Column Title:',font=font_type,bg='#f5f5f5').place(x=45,y=70)
            Label(top,text='Row Identifier:',font=font_type,bg='#f5f5f5').place(x=45,y=130)
            Label(top,text='Values:',font=font_type,bg='#f5f5f5').place(x=45,y=190)
            Label(top,text='Data aggregation:',font=font_type,bg='#f5f5f5').place(x=10,y=250)
            aggregation=Entry(top,width=15,font=font_type,bg='#fcfcfc',relief='groove')
            aggregation.place(x=115,y=250)
            aggregation.insert(INSERT,'median')
            aggregation.config(state=DISABLED,disabledbackground='#fcfcfc')
            Button(top,image=foto25,bg='#f5f5f5',relief='flat',overrelief='groove',cursor='hand2',command=aggregate_method).place(x=210,y=245)
            Button(top,text='Help',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2').place(x=10,y=310)
            Button(top,text='OK',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=start_pivot).place(x=200,y=310)
            Button(top,text='Cancel',width=7,font=font_type,bg='#fcfcfc',relief='groove',cursor='hand2',command=cancel_pivot).place(x=265,y=310)
        if len(file_name)>=1:
            warning=showinfo(title='Data Pivoting',message='Save and close this project before opening a new one')
    #Descripcion de todos los elementos presentes en el menu de la aplicacion
    menubar=Menu(win)
    filemenu=Menu(menubar,tearoff=0,font=font_type,bg='#fcfcfc')
    editmenu=Menu(menubar,tearoff=0,font=font_type,bg='#fcfcfc')
    toolsmenu=Menu(menubar,tearoff=0,font=font_type,bg='#fcfcfc')
    technologies=Menu(menubar,tearoff=0,font=font_type,bg='#fcfcfc')
    helpmenu=Menu(menubar,tearoff=0,font=font_type,bg='#fcfcfc')
    microarrays=Menu(filemenu,tearoff=0,font=font_type,bg='#fcfcfc')
    one_channel=Menu(microarrays,tearoff=0,font=font_type,bg='#fcfcfc')
    two_channel=Menu(microarrays,tearoff=0,font=font_type,bg='#fcfcfc')
    openmenu=Menu(microarrays,tearoff=0,font=font_type,bg='#fcfcfc')
    HTA=Menu(filemenu,tearoff=0,font=font_type,bg='#fcfcfc')
    preprocessing=Menu(toolsmenu,tearoff=0,font=font_type,bg='#fcfcfc')
    statistics=Menu(toolsmenu,tearoff=0,font=font_type,bg='#fcfcfc')
    functional=Menu(toolsmenu,tearoff=0,font=font_type,bg='#fcfcfc')
    one_channel.add_command(label='Affymetrix Microarrays',command=affymetrix_technology)
    one_channel.add_command(label='Illumina Microarrays',command=illumina_technology)
    one_channel.add_command(label='Custom Microarrays',command=custom_technology)
    two_channel.add_command(label='Agilent Microarrays',command=agilent_technology)
    two_channel.add_command(label='Genepix Microarrays',command=genepix_technology)
    two_channel.add_command(label='Imagene Microarrays',command=imagene_technology)
    two_channel.add_command(label='Quantarray Microarrays',command=quantarray_technology)
    two_channel.add_command(label='Scanarray Microarrays',command=scanarray_technology)
    preprocessing.add_command(label='Normalization',command=normalize_data)
    preprocessing.add_command(label='Imputation',command=impute_data)
    preprocessing.add_command(label='Summarize Data',command=summarize_data)
    preprocessing.add_command(label='Transpose Data',command=transpose_data)
    statistics.add_command(label='Linear Methods for Microarrays (LIMMA)',command=limma_data_analysis)
    statistics.add_command(label='Significant Analisys of Microarrays (SAM)',command=sam_data_analysis)
    statistics.add_command(label='Unpaired T-Test',command=ttest_data_analysis)
    functional.add_command(label='Hierarchical Clustering',command=cluster_data)
    functional.add_command(label='Retrieve Annotations',command=annotation_retriever)
    filemenu.add_cascade(label='Open',menu=openmenu)
    filemenu.add_separator()
    filemenu.add_command(label='Save',command=save_file)
    filemenu.add_command(label='Save as',command=export_project)
    filemenu.add_separator()
    filemenu.add_command(label='Clear',command=start_new_project)
    filemenu.add_separator()
    filemenu.add_command(label='Exit',command=exit_genomics_office)
    openmenu.add_command(label='Open',command=open_file)
    openmenu.add_command(label='Open GQ project', command=open_GeneQuery_project)
    openmenu.add_command(label='Open + Pivot',command=pivot_data)
    editmenu.add_command(label='Copy (Ctrl + C)',command=clipboard)
    editmenu.add_command(label='Paste (Ctrl + V)',command=paste_clipboard)
    toolsmenu.add_cascade(label='Preprocessing',menu=preprocessing)
    toolsmenu.add_cascade(label='Statistical Analysis',menu=statistics)
    toolsmenu.add_cascade(label='Clustering and Functional analysis',menu=functional)
    microarrays.add_cascade(label='One channel technologies',menu=one_channel)
    microarrays.add_cascade(label='Two channel technologies',menu=two_channel)
    HTA.add_command(label='qPCR',command=qPCR_technology)
    HTA.add_command(label='CGH analysis',command=CGH_technology)
    technologies.add_cascade(label='Microarrays',menu=microarrays)
    technologies.add_cascade(label='High Throughput Assays',menu=HTA)
    helpmenu.add_command(label='PDF Manual',command=manuals)
    helpmenu.add_command(label='Getting Started',command=user_case)
    helpmenu.add_separator()
    helpmenu.add_command(label='GeneQuery Support',command=about_GeneQuery)
    helpmenu.add_separator()
    helpmenu.add_command(label='About GeneQuery',command=about_GeneQuery)
    menubar.add_cascade(label='File',menu=filemenu)
    menubar.add_cascade(label='Edit',menu=editmenu)
    menubar.add_cascade(label='Tools',menu=toolsmenu)
    menubar.add_cascade(label='Workflows',menu=technologies)
    menubar.add_cascade(label='Help',menu=helpmenu)
    win.config(menu=menubar)
    win.protocol("WM_DELETE_WINDOW",exit_genomics_office)
    win.bind("<Control-c>",clipboard_callback)
    win.bind("<Control-v>",paste_callback)
    win.mainloop()
#Imagenes de los botones situados en la barra de herramientas
foto1=PhotoImage(file='C://PG//Icons//Folder.gif')
foto2=PhotoImage(file='C://PG//Icons//Save.gif')
foto3=PhotoImage(file='C://PG//Icons//Printer.gif')
foto4=PhotoImage(file='C://PG//Icons//Document.gif')
foto5=PhotoImage(file='C://PG//Icons//Clipboard Copy.gif')
foto6=PhotoImage(file='C://PG//Icons//Clipboard Cut.gif')
foto7=PhotoImage(file='C://PG//Icons//Clipboard Paste.gif')
foto8=PhotoImage(file='C://PG//Icons//Dots.gif')
foto9=PhotoImage(file='C://PG//Icons//Tool.gif')
foto10=PhotoImage(file='C://PG//Icons//Go In.gif')
foto11=PhotoImage(file='C://PG//Icons//Light.gif')
foto12=PhotoImage(file='C://PG//Icons//Stats.gif')
foto13=PhotoImage(file='C://PG//Icons//Poll.gif')
foto14=PhotoImage(file='C://PG//Icons//Direction Horz.gif')
foto15=PhotoImage(file='C://PG//Icons//Document New.gif')
foto16=PhotoImage(file='C://PG//Icons//Question.gif')
foto17=PhotoImage(file='C://PG//Icons//Tick.gif')
foto18=PhotoImage(file='C://PG//Icons//Wrong.gif')
foto19=PhotoImage(file='C://PG//Icons//Zoom In.gif')
foto20=PhotoImage(file='C://PG//Icons//Zoom Out.gif')
foto21=PhotoImage(file='C://PG//Icons//Gear.gif')
foto22=PhotoImage(file='C://PG//Icons//Go Out.gif')
foto23=PhotoImage(file='C://PG//Icons//Refresh.gif')
foto24=PhotoImage(file='C://PG//Icons//Arrow1 Up.gif')
foto25=PhotoImage(file='C://PG//Icons//Arrow1 Down.gif')
foto26=PhotoImage(file='C://PG//Icons//Arrow1 Right.gif')
foto27=PhotoImage(file='C://PG//Icons//Arrow1 Left.gif')
foto28=PhotoImage(file='C://PG//Icons//Logo.gif')
foto29=PhotoImage(file='C://PG//Icons//Application.gif')
foto30=PhotoImage(file='C://PG//Icons//Applications.gif')
foto31=PhotoImage(file='C://PG//Icons//Dots Down.gif')
foto32=PhotoImage(file='C://PG//Icons//Box.gif')
foto33=PhotoImage(file='C://PG//Icons//Directions.gif')
foto34=PhotoImage(file='C://PG//Icons//Sitemap.gif')
foto35=PhotoImage(file='C://PG//Icons//Write2.gif')
foto36=PhotoImage(file='C://PG//Icons//Pivot1.gif')
foto37=PhotoImage(file='C://PG//Icons//Pivot2.gif')
foto38=PhotoImage(file='C://PG//Icons//Pivot3.gif')
new_project()
