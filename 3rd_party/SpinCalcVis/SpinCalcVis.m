function [] = SpinCalcVis()
%==========================================================================
%SpinCalcVis Rotational Visualization Tool  V1.0
%Author: John Fuller
%        imaginationdirect@gmail.com
%
%This is an instructional GUI to be used for learning how Euler angles,
%DCMs, quaternions, and Euler vector parameters relate in rotation of
%cartesian frames (A to B).  
%
%For the function-based rotation conversion, please see SpinCalc:
%
%http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-conv
%ert-between-dcm-euler-angles-quaternions-and-euler-vectors
%
%Uses an enhanced uicontrol GUI function for support of LaTeX formatting:
%Function uibutton, Author: Douglas Schwarz
%
%Source:
%http://www.mathworks.com/matlabcentral/fileexchange/10743-uibutton-gui-pus
%hbuttons-with-better-labels
%==========================================================================
S.fh = figure('units','pixels',...
              'position',[300 150 940 500],...
              'menubar','none',...
              'resize','off',...
              'numbertitle','off',...
              'name','SpinCalcVis: Rotation Visualization GUI');

S.axes1=axes('units','pixels',...
             'position',[460 12 470 470],...
             'view',[58 28],...
             'xtick',[],...
             'ytick',[],...
             'ztick',[]);
         
plot3([0;1],[0;0],[0;0],':r')
hold on
axis equal
text(1.05,0,-0.01,'X_{A}','color','r')
plot3([0;0],[0;1],[0;0],':g')
text(0,1.05,-0.01,'Y_{A}','color','g')
plot3([0;0],[0;0],[0;1],':b')
text(-0.01,-0.01,1.05,'Z_{A}','color','b')
xlim([-1.0,1.0])
ylim([-1.0,1.0])
zlim([-1.0,1.0])
set(gca,'xtick',[],'ytick',[],'ztick',[],'box','off')
view([58 28])
hold off
rotate3d(S.axes1)          

S.fr1 = uicontrol('style','frame',...
                 'units','pixels',...
                 'position',[10 405 440 85]);
            
S.title1 = uicontrol('style','text',...
                    'units','pixels',...
                    'position',[(450)/2-100 460-20 200 40],...
                    'string',{'SpinCalcVis';'Rotational Visualization Tool'},...
                    'FontSize',10,...
                    'FontWeight','bold');
                
S.title2 = uicontrol('style','text',...
                    'units','pixels',...
                    'position',get(S.title1,'position')+[0 -34 0 -5],...
                    'string',{'Author: John Fuller';'Version 1.0'});
                 
S.fr2 = uicontrol('style','frame',...
                 'units','pixels',...
                 'position',[20 10 420 60]);

S.fr3 = uicontrol('style','frame',...
                 'units','pixels',...
                 'position',[358 255 90 58]);             

S.txtradios = uicontrol('style','text',...
                        'units','pixels',...
                        'position',[358 320 80 12],...
                        'BackgroundColor',0.7967*[1 1 1],...
                        'string','Plotting');             
                    
S.dialog = uicontrol('style','text',...
                    'units','pixels',...
                    'position',[24 12 415 56],...
                    'HorizontalAlignment','left',...
                    'Foregroundcolor',0.3*[1 1 1],...
                    'string',{'No Output'});
                
S.dialogtxt = uicontrol('style','text',...
                     'units','pixels',...
                     'position',[12, 70 60 18],...
                     'BackgroundColor',0.7967*[1 1 1],...
                     'string','Dialog Box');
                
S.txt1 = uicontrol('style','text',...
                   'units','pixels',...
                   'position',[13 380 100 12],...
				   'BackgroundColor',0.7967*[1 1 1],...
                   'string','Euler Angles (deg)');
               
S.ea1 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 355 80 18],...
                  'BackgroundColor',[1 1 1],...
                  'String',0,...
                  'Value',0);

S.ea2 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.ea1,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);

S.ea3 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.ea2,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);

S.ea1txt = uibutton('style','text',...
                    'units','pixels',...
                    'position',[20 338 80 18],...
				    'BackgroundColor',0.7967*[1 1 1],...
                    'string','\psi');              

S.ea2txt = uibutton('style','text',...
                    'units','pixels',...
                    'position',[110 338 80 18],...
				    'BackgroundColor',0.7967*[1 1 1],...
                    'string','\theta');
                
S.ea3txt = uibutton('style','text',...
                    'units','pixels',...
                    'position',[200 338 80 18],...
				    'BackgroundColor',0.7967*[1 1 1],...
                    'string','\phi');  

set(S.ea1,'Callback',{@assignea1,S});
set(S.ea2,'Callback',{@assignea2,S});
set(S.ea3,'Callback',{@assignea3,S});              
              
S.eamenu = uicontrol('style','popupmenu',...
                     'units','pixels',...
                     'position',get(S.ea3,'position')+[90 2 0 0],...
                     'BackgroundColor',[1 1 1],...
                     'String',{'XYZ  (123)';'XZY  (132)';'YXZ  (213)';'YZX  (231)';'ZXY  (312)';'ZYX  (321)';'XYX  (121)';'XZX  (131)';'YXY  (212)';'YZY  (232)';'ZXZ  (313)';'ZYZ  (323)'});
                 
S.txt2 = uicontrol('style','text',...
                   'units','pixels',...
                   'position',get(S.eamenu,'position')+[-2 18 0 -5],...
				   'BackgroundColor',0.7967*[1 1 1],...
                   'string','Order');
			   
S.txt3 = uicontrol('style','text',...
                   'units','pixels',...
                   'position',[13 320 160 12],...
				   'BackgroundColor',0.7967*[1 1 1],...
                   'string','Directions Cosine Matrix (DCM)');
			   
S.dcm11 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 295 80 18],...
                  'String',1,...
                  'Value',1,...
                  'BackgroundColor',[1 1 1]);
			 
S.dcm12 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm11,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.dcm13 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm12,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.dcm21 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 275 80 18],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			 
S.dcm22 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm21,'position')+[90 0 0 0],...
                  'String',1,...
                  'Value',1,...
                  'BackgroundColor',[1 1 1]);
			  
S.dcm23 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm22,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.dcm31 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 255 80 18],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			 
S.dcm32 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm31,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.dcm33 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.dcm32,'position')+[90 0 0 0],...
                  'String',1,...
                  'Value',1,...
                  'BackgroundColor',[1 1 1]);

set(S.dcm11,'Callback',{@assigndcm11,S});
set(S.dcm12,'Callback',{@assigndcm12,S});
set(S.dcm13,'Callback',{@assigndcm13,S});
set(S.dcm21,'Callback',{@assigndcm21,S});
set(S.dcm22,'Callback',{@assigndcm22,S});
set(S.dcm23,'Callback',{@assigndcm23,S});
set(S.dcm31,'Callback',{@assigndcm31,S});
set(S.dcm32,'Callback',{@assigndcm32,S});
set(S.dcm33,'Callback',{@assigndcm33,S});
              
S.txt4 = uicontrol('style','text',...
                   'units','pixels',...
                   'position',[13 225 64 12],...
				   'BackgroundColor',0.7967*[1 1 1],...
                   'string','Quaternion');
			   
S.q1 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 200 80 18],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
              
S.q2 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.q1,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);

S.q3 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.q2,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.q4 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.q3,'position')+[90 0 0 0],...
                  'String',1,...
                  'Value',1,...
                  'BackgroundColor',[1 1 1]);

set(S.q1,'Callback',{@assignq1,S});
set(S.q2,'Callback',{@assignq2,S});
set(S.q3,'Callback',{@assignq3,S});
set(S.q4,'Callback',{@assignq4,S});              
              
S.q1txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[20 180 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{1}sin(\mu/2)');

S.q2txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[110 180 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{2}sin(\mu/2)');
			  
S.q3txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[200 180 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{3}sin(\mu/2)');
			  
S.q3txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[290 183 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','cos(\mu/2)');
			  
S.txt6 = uicontrol('style','text',...
                   'units','pixels',...
                   'position',[13 155 90 12],...
				   'BackgroundColor',0.7967*[1 1 1],...
                   'string','Euler Parameters');
			   
S.ep1 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',[20 130 80 18],...
                  'String',1,...
                  'Value',1,...
                  'BackgroundColor',[1 1 1]);             
              
S.ep2 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.ep1,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);

S.ep3 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.ep2,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);
			  
S.ep4 = uicontrol('style','edit',...
                  'units','pixels',...
                  'position',get(S.ep3,'position')+[90 0 0 0],...
                  'String',0,...
                  'Value',0,...
                  'BackgroundColor',[1 1 1]);

set(S.ep1,'Callback',{@assignep1,S});
set(S.ep2,'Callback',{@assignep2,S});
set(S.ep3,'Callback',{@assignep3,S});
set(S.ep4,'Callback',{@assignep4,S});
              
S.ep1txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[20 110 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{1}');

S.ep2txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[110 110 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{2}');
			  
S.ep3txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[200 110 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','m_{3}');
			  
S.ep3txt = uibutton('style','text',...
                  'units','pixels',...
                  'position',[290 115 80 18],...
				  'BackgroundColor',0.7967*[1 1 1],...
                  'string','\mu');

S.exitbutton = uibutton('style','pushbutton',...
                        'units','pixels',...
                        'position',[320 80 120 25],...
                        'callback',{@closegui,S},...
                        'string','Exit');


S.plotea = uicontrol('style','radio',...
                     'units','pixels',...
                     'position',[362 290 80 14],...
                     'string','Euler Angles',...
                     'fontsize',8,...
                     'value',1);

        
                 
S.plotep = uicontrol('style','radio',...
                     'units','pixels',...
                     'position',[362 268 80 14],...
                     'string','Euler Vector',...
                     'fontsize',8,...
                     'value',1); 

set(S.plotea,'callback',{@replot,S});                   
set(S.plotep,'callback',{@replot,S}); 
                    
S.pbea = uicontrol('style','pushbutton',...
                   'units','pixels',...
                   'position',get(S.eamenu,'position')+[90 -2 -20 0],...
                   'callback',{@spincalcea,S},...
                   'String','Convert');
              

S.pbq = uicontrol('style','pushbutton',...
                   'units','pixels',...
                   'position',get(S.q4,'position')+[90 0 -20 0],...
                   'callback',{@spincalcq,S},...
                   'String','Convert');
               
S.pbnormq = uicontrol('style','pushbutton',...
                      'units','pixels',...
                      'position',get(S.pbq,'position')+[0 20 0 0],...
                      'callback',{@normq,S},...
                      'String','Normalize');                  
                  
S.pbep = uicontrol('style','pushbutton',...
                   'units','pixels',...
                   'position',get(S.ep4,'position')+[90 0 -20 0],...
                   'callback',{@spincalcep,S},...
                   'String','Convert');
          
S.pbnormep = uicontrol('style','pushbutton',...
                       'units','pixels',...
                       'position',get(S.pbep,'position')+[0 20 0 0],...
                       'callback',{@normep,S},...
                       'String','Normalize');
              
S.pbdcm = uicontrol('style','pushbutton',...
                    'units','pixels',...
                    'position',get(S.dcm33,'position')+[90 0 -20 0],...
                    'callback',{@spincalcdcm,S},...
                    'String','Convert');
               
S.pbtransdcm = uicontrol('style','pushbutton',...
                         'units','pixels',...
                         'position',get(S.pbdcm,'position')+[0 20 0 0],...
                         'callback',{@transdcm,S},...
                         'String','Transpose');                   
end

function []=replot(varargin)
    S=varargin{3};
    S=plotter(S);
end

function S=plotter(S)
    axes(S.axes1);
    currentview=get(S.axes1,'view');
    plot3([0;1],[0;0],[0;0],':r')
    hold on
    axis equal
    text(1.05,0,-0.01,'X_{A}','color','r')
    plot3([0;0],[0;1],[0;0],':g')
    text(0,1.05,-0.01,'Y_{A}','color','g')
    plot3([0;0],[0;0],[0;1],':b')
    text(-0.01,-0.01,1.05,'Z_{A}','color','b')
    xlim([-1.0,1.0])
    ylim([-1.0,1.0])
    zlim([-1.0,1.0])
    set(gca,'xtick',[],'ytick',[],'ztick',[],'box','off')    
    xlim([-1.0,1.0])
    ylim([-1.0,1.0])
    zlim([-1.0,1.0])
    dcm11=get(S.dcm11,'value');
    dcm12=get(S.dcm12,'value');
    dcm13=get(S.dcm13,'value');
    dcm21=get(S.dcm21,'value');
    dcm22=get(S.dcm22,'value');
    dcm23=get(S.dcm23,'value');
    dcm31=get(S.dcm31,'value');
    dcm32=get(S.dcm32,'value');
    dcm33=get(S.dcm33,'value');
    S.XB=plot3([0;dcm11],[0;dcm12],[0;dcm13],'-r');
    S.XBtext=text(dcm11*1.05,dcm12*1.05,dcm13*1.05,'X_{B}','color','r');
    S.YB=plot3([0;dcm21],[0;dcm22],[0;dcm23],'-g');
    S.YBtext=text(dcm21*1.05,dcm22*1.05,dcm23*1.05,'Y_{B}','color','g');
    S.ZB=plot3([0;dcm31],[0;dcm32],[0;dcm33],'-b');
    S.ZBtext=text(dcm31*1.05,dcm32*1.05,dcm33*1.05,'Z_{B}','color','b');
    if get(S.plotea,'value')==1               
        ea1=get(S.ea1,'value');
        ea2=get(S.ea2,'value');
        ea3=get(S.ea3,'value');
        EA_rotation_order_index=get(S.eamenu,'value');
        %EA_rotations={'XYZ';'XZY';'YXZ';'YZX';'ZXY';'ZYX';'XYX';'XZX';'YXY';'YZY';'ZXZ';'ZYZ'};
        EA_rotations={'EA123';'EA132';'EA213';'EA231';'EA312';'EA321';'EA121';'EA131';'EA212';'EA232';'EA313';'EA323'};
        EA_rotation_order=EA_rotations{EA_rotation_order_index,:};

        n1=max([ceil(abs(ea1)/0.5),19]);
        n2=max([ceil(abs(ea2)/0.5),19]);
        n3=max([ceil(abs(ea3)/0.5),19]);
        EA_sweep1=[linspace(0,ea1,n1)',zeros(n1,1),zeros(n1,1)];
        EA_sweep2=[ea1*ones(n2,1),linspace(0,ea2,n2)',zeros(n2,1)];
        EA_sweep3=[ea1*ones(n3,1),ea2*ones(n3,1),linspace(0,ea3,n3)'];
        full_sweep=[EA_sweep1;EA_sweep2;EA_sweep3];
        if EA_rotation_order_index<=6
            [DCM_sweep,errorstring]=spincalcmod([EA_rotation_order,'toDCM'],full_sweep,eps,0);
        else
            EA_rotation_order_mod=EA_rotations{EA_rotation_order_index-6,:};
            [DCM_sweep1,errorstring]=spincalcmod([EA_rotation_order_mod,'toDCM'],full_sweep(1:(n1+1),1:3),eps,0);
            [DCM_sweep2,errorstring]=spincalcmod([EA_rotation_order,'toDCM'],full_sweep((n1+2):end,1:3),eps,0);
            DCM_sweep=cat(3,DCM_sweep1,DCM_sweep2);
        end
        sweep_vector1a=NaN(n1,3);
        sweep_vector1b=NaN(n1,3);
        sweep_vector2a=NaN(n2,3);
        sweep_vector2b=NaN(n2,3);
        sweep_vector3a=NaN(n3,3);
        sweep_vector3b=NaN(n3,3);
        initvec1a=[0,0,0];
        initvec1b=[0,0,0];
        initvec2a=[0,0,0];
        initvec2b=[0,0,0];
        initvec3a=[0,0,0];
        initvec3b=[0,0,0];
        colors={'r';'g';'b'};
        angle1_color=colors{str2double(EA_rotations{EA_rotation_order_index}(3))};
        angle2_color=colors{str2double(EA_rotations{EA_rotation_order_index}(4))};
        angle3_color=colors{str2double(EA_rotations{EA_rotation_order_index}(5))};
        if EA_rotation_order_index<=6
            eval(['initvec1a(1,',EA_rotations{EA_rotation_order_index}(4),')=0.8;colorsweep1a=colors{',EA_rotations{EA_rotation_order_index}(4),'};'])
            eval(['initvec1b(1,',EA_rotations{EA_rotation_order_index}(5),')=0.5;colorsweep1b=colors{',EA_rotations{EA_rotation_order_index}(5),'};'])
            eval(['initvec2a(1,',EA_rotations{EA_rotation_order_index}(3),')=0.65;colorsweep2a=colors{',EA_rotations{EA_rotation_order_index}(3),'};'])
            eval(['initvec2b(1,',EA_rotations{EA_rotation_order_index}(5),')=0.5;'])
            eval(['initvec3a(1,',EA_rotations{EA_rotation_order_index}(3),')=0.65;'])
            eval(['initvec3b(1,',EA_rotations{EA_rotation_order_index}(4),')=0.8;'])
        else
            if str2double(EA_rotations{EA_rotation_order_index}(3))==1
                if str2double(EA_rotations{EA_rotation_order_index}(4))==2
                    other_axis=3; %#ok<*NASGU>
                else
                    other_axis=2;
                end
            elseif str2double(EA_rotations{EA_rotation_order_index}(3))==2
                if str2double(EA_rotations{EA_rotation_order_index}(4))==1
                    other_axis=3;
                else
                    other_axis=1;
                end
            elseif str2double(EA_rotations{EA_rotation_order_index}(3))==3
                if str2double(EA_rotations{EA_rotation_order_index}(4))==1
                    other_axis=2;
                else
                    other_axis=1;
                end
            end
            eval(['initvec1a(1,',EA_rotations{EA_rotation_order_index}(4),')=0.8;colorsweep1a=colors{',EA_rotations{EA_rotation_order_index}(4),'};'])
            eval('initvec1b(1,other_axis)=0.5;colorsweep1b=colors{other_axis};')
            eval(['initvec2a(1,',EA_rotations{EA_rotation_order_index}(3),')=0.65;colorsweep2a=colors{',EA_rotations{EA_rotation_order_index}(3),'};'])
            eval('initvec2b(1,other_axis)=0.5;')
            eval('initvec3a(1,other_axis)=0.5;')
            eval(['initvec3b(1,',EA_rotations{EA_rotation_order_index}(4),')=0.8;'])
        end
        for ii=1:n1
            sweep_vector1a(ii,1:3)=initvec1a*DCM_sweep(1:3,1:3,ii);
            sweep_vector1b(ii,1:3)=initvec1b*DCM_sweep(1:3,1:3,ii);
        end
        for ii=1:n2
            sweep_vector2a(ii,1:3)=initvec2a*DCM_sweep(1:3,1:3,n1+ii);
            sweep_vector2b(ii,1:3)=initvec2b*DCM_sweep(1:3,1:3,n1+ii);
        end
        for ii=1:n3
            sweep_vector3a(ii,1:3)=initvec3a*DCM_sweep(1:3,1:3,n1+n2+ii);
            sweep_vector3b(ii,1:3)=initvec3b*DCM_sweep(1:3,1:3,n1+n2+ii);
        end
        if abs(ea1)<1e-10
            sweep_vector1a(:,:)=NaN;
            sweep_vector1b(:,:)=NaN;
        end
        if abs(ea2)<1e-10
            sweep_vector2a(:,:)=NaN;
            sweep_vector2b(:,:)=NaN;
        end
        if abs(ea3)<1e-10
            sweep_vector3a(:,:)=NaN;
            sweep_vector3b(:,:)=NaN;
        end
        %Determine coordinates of arrowtips
        for ii=1:6
            if ii==1,
                temp=sweep_vector1a;
            elseif ii==2
                temp=sweep_vector1b;
            elseif ii==3
                temp=sweep_vector2a;
            elseif ii==4
                temp=sweep_vector2b;
            elseif ii==5
                temp=sweep_vector3a;
            elseif ii==6
                temp=sweep_vector3b;
            end
            Pend=temp(end,1:3);
            Pstart=temp(end-18,1:3);
            Pmid=temp(end-9,1:3);
            Pnear=(Pstart+Pend)/2;
            Pa=8*(Pnear-Pmid)+Pmid;
            Pb=-8*(Pnear-Pmid)+Pmid;
            x=[Pend(1);Pa(1);Pb(1);Pend(1)];
            y=[Pend(2);Pa(2);Pb(2);Pend(2)];
            z=[Pend(3);Pa(3);Pb(3);Pend(3)];
            S.arrowheads=plot3(x,y,z,'-k');
        end
        S.start1a=plot3(sweep_vector1a(1,1),sweep_vector1a(1,2),sweep_vector1a(1,3),'ok','markersize',5);
        S.sweep1a=plot3(sweep_vector1a(1:end-9,1),sweep_vector1a(1:end-9,2),sweep_vector1a(1:end-9,3),'-k');        
        S.sweep1atick=plot3([sweep_vector1a(ceil(n1/2),1);sweep_vector1a(ceil(n1/2),1)/0.92],[sweep_vector1a(ceil(n1/2),2);sweep_vector1a(ceil(n1/2),2)/0.92],[sweep_vector1a(ceil(n1/2),3);sweep_vector1a(ceil(n1/2),3)/0.92],'-k');
        S.psi_a=text(sweep_vector1a(ceil(n1/2),1)/0.9,sweep_vector1a(ceil(n1/2),2)/0.9,sweep_vector1a(ceil(n1/2),3)/0.9,'\psi','fontsize',11,'color',angle1_color);
        S.axis1a=plot3([0,sweep_vector1a(end,1)],[0,sweep_vector1a(end,2)],[0,sweep_vector1a(end,3)],[':',colorsweep1a]);
        
        S.start1b=plot3(sweep_vector1b(1,1),sweep_vector1b(1,2),sweep_vector1b(1,3),'ok','markersize',5);
        S.sweep1b=plot3(sweep_vector1b(1:end-9,1),sweep_vector1b(1:end-9,2),sweep_vector1b(1:end-9,3),'-k');        
        S.sweep1btick=plot3([sweep_vector1b(ceil(n1/2),1);sweep_vector1b(ceil(n1/2),1)/0.92],[sweep_vector1b(ceil(n1/2),2);sweep_vector1b(ceil(n1/2),2)/0.92],[sweep_vector1b(ceil(n1/2),3);sweep_vector1b(ceil(n1/2),3)/0.92],'-k');
        S.psi_b=text(sweep_vector1b(ceil(n1/2),1)/0.9,sweep_vector1b(ceil(n1/2),2)/0.9,sweep_vector1b(ceil(n1/2),3)/0.9,'\psi','fontsize',11,'color',angle1_color);
        S.axis1b=plot3([0,sweep_vector1b(end,1)],[0,sweep_vector1b(end,2)],[0,sweep_vector1b(end,3)],[':',colorsweep1b]);
        
        S.start2a=plot3(sweep_vector2a(1,1),sweep_vector2a(1,2),sweep_vector2a(1,3),'ok','markersize',5);
        S.sweep2a=plot3(sweep_vector2a(1:end-9,1),sweep_vector2a(1:end-9,2),sweep_vector2a(1:end-9,3),'-k');
        S.sweep2atick=plot3([sweep_vector2a(ceil(n2/2),1);sweep_vector2a(ceil(n2/2),1)/0.92],[sweep_vector2a(ceil(n2/2),2);sweep_vector2a(ceil(n2/2),2)/0.92],[sweep_vector2a(ceil(n2/2),3);sweep_vector2a(ceil(n2/2),3)/0.92],'-k');
        S.theta_a=text(sweep_vector2a(ceil(n2/2),1)/0.9,sweep_vector2a(ceil(n2/2),2)/0.9,sweep_vector2a(ceil(n2/2),3)/0.9,'\theta','fontsize',11,'color',angle2_color);
        S.axis2a=plot3([0,sweep_vector2a(end,1)],[0,sweep_vector2a(end,2)],[0,sweep_vector2a(end,3)],[':',colorsweep2a]);
        
        S.start2b=plot3(sweep_vector2b(1,1),sweep_vector2b(1,2),sweep_vector2b(1,3),'ok','markersize',5);
        S.sweep2b=plot3(sweep_vector2b(1:end-9,1),sweep_vector2b(1:end-9,2),sweep_vector2b(1:end-9,3),'-k');
        S.sweep2btick=plot3([sweep_vector2b(ceil(n2/2),1);sweep_vector2b(ceil(n2/2),1)/0.92],[sweep_vector2b(ceil(n2/2),2);sweep_vector2b(ceil(n2/2),2)/0.92],[sweep_vector2b(ceil(n2/2),3);sweep_vector2b(ceil(n2/2),3)/0.92],'-k');
        S.theta_b=text(sweep_vector2b(ceil(n2/2),1)/0.9,sweep_vector2b(ceil(n2/2),2)/0.9,sweep_vector2b(ceil(n2/2),3)/0.9,'\theta','fontsize',11,'color',angle2_color);
        
        S.start3a=plot3(sweep_vector3a(1,1),sweep_vector3a(1,2),sweep_vector3a(1,3),'ok','markersize',5);
        S.sweep3a=plot3(sweep_vector3a(1:end-9,1),sweep_vector3a(1:end-9,2),sweep_vector3a(1:end-9,3),'-k');
        S.sweep3atick=plot3([sweep_vector3a(ceil(n3/2),1);sweep_vector3a(ceil(n3/2),1)/0.92],[sweep_vector3a(ceil(n3/2),2);sweep_vector3a(ceil(n3/2),2)/0.92],[sweep_vector3a(ceil(n3/2),3);sweep_vector3a(ceil(n3/2),3)/0.92],'-k');
        S.phi_a=text(sweep_vector3a(ceil(n3/2),1)/0.9,sweep_vector3a(ceil(n3/2),2)/0.9,sweep_vector3a(ceil(n3/2),3)/0.9,'\phi','fontsize',11,'color',angle3_color);
        if EA_rotation_order_index>6
            S.axis3a=plot3([0,sweep_vector2b(end,1)],[0,sweep_vector2b(end,2)],[0,sweep_vector2b(end,3)],[':',colorsweep1b]);
        end
        
        S.start3b=plot3(sweep_vector3b(1,1),sweep_vector3b(1,2),sweep_vector3b(1,3),'ok','markersize',5);
        S.sweep3b=plot3(sweep_vector3b(1:end-9,1),sweep_vector3b(1:end-9,2),sweep_vector3b(1:end-9,3),'-k');
        S.sweep3btick=plot3([sweep_vector3b(ceil(n3/2),1);sweep_vector3b(ceil(n3/2),1)/0.92],[sweep_vector3b(ceil(n3/2),2);sweep_vector3b(ceil(n3/2),2)/0.92],[sweep_vector3b(ceil(n3/2),3);sweep_vector3b(ceil(n3/2),3)/0.92],'-k');
        S.phi_b=text(sweep_vector3b(ceil(n3/2),1)/0.9,sweep_vector3b(ceil(n3/2),2)/0.9,sweep_vector3b(ceil(n3/2),3)/0.9,'\phi','fontsize',11,'color',angle3_color);
    end
    if get(S.plotep,'value')==1
        m1=get(S.ep1,'value');
        m2=get(S.ep2,'value');
        m3=get(S.ep3,'value');
        plot3([0 m1],[0 m2],[0 m3],'-m');
        text(m1*1.02,m2*1.02,m3*1.02,'m','color','m');
    end
    hold off
    set(gca,'xtick',[],'ytick',[],'ztick',[],'box','off')
    view(currentview)
end

function []=normq(varargin)
    S=varargin{3};
    q1=get(S.q1,'value');
	q2=get(S.q2,'value');
	q3=get(S.q3,'value');
	q4=get(S.q4,'value');
    qnorm=norm([q1,q2,q3,q4]);
    set(S.q1,'string',num2str(q1/qnorm,'%8.6f'),'value',q1/qnorm);
    set(S.q2,'string',num2str(q2/qnorm,'%8.6f'),'value',q2/qnorm);
    set(S.q3,'string',num2str(q3/qnorm,'%8.6f'),'value',q3/qnorm);
    set(S.q4,'string',num2str(q4/qnorm,'%8.6f'),'value',q4/qnorm);
end

function []=normep(varargin)
    S=varargin{3};
    m1=get(S.ep1,'value');
	m2=get(S.ep2,'value');
	m3=get(S.ep3,'value');
    mnorm=norm([m1,m2,m3]);
    set(S.ep1,'string',num2str(m1/mnorm,'%8.6f'),'value',m1/mnorm);
    set(S.ep2,'string',num2str(m2/mnorm,'%8.6f'),'value',m2/mnorm);
    set(S.ep3,'string',num2str(m3/mnorm,'%8.6f'),'value',m3/mnorm);
end

function []=transdcm(varargin)
    S=varargin{3};
	dcm12=get(S.dcm12,'value');
	dcm13=get(S.dcm13,'value');
    dcm21=get(S.dcm21,'value');
	dcm23=get(S.dcm23,'value');
    dcm31=get(S.dcm31,'value');
	dcm32=get(S.dcm32,'value');
    set(S.dcm12,'string',num2str(dcm21,'%8.6f'),'value',dcm21);
    set(S.dcm13,'string',num2str(dcm31,'%8.6f'),'value',dcm31);
    set(S.dcm21,'string',num2str(dcm12,'%8.6f'),'value',dcm12);
    set(S.dcm23,'string',num2str(dcm32,'%8.6f'),'value',dcm32);
    set(S.dcm31,'string',num2str(dcm13,'%8.6f'),'value',dcm13);
    set(S.dcm32,'string',num2str(dcm23,'%8.6f'),'value',dcm23);
end

function []=spincalcea(varargin)
	S=varargin{3};
	ea1=get(S.ea1,'value');
	ea2=get(S.ea2,'value');
	ea3=get(S.ea3,'value');
	EA_rotation_order_index=get(S.eamenu,'value');
	EA_rotations={'EA123';'EA132';'EA213';'EA231';'EA312';'EA321';'EA121';'EA131';'EA212';'EA232';'EA313';'EA323'};
	EA_rotation_order=EA_rotations{EA_rotation_order_index,:};
	[DCM,errorstring]=spincalcmod([EA_rotation_order,'toDCM'],[ea1,ea2,ea3],eps,1);
	[EV,errorstring]=spincalcmod([EA_rotation_order,'toEV'],[ea1,ea2,ea3],eps,1);
	[Q,errorstring]=spincalcmod([EA_rotation_order,'toQ'],[ea1,ea2,ea3],eps,1);
	set(S.dcm11,'string',num2str(DCM(1,1),'%8.6f'),'value',DCM(1,1));
	set(S.dcm12,'string',num2str(DCM(1,2),'%8.6f'),'value',DCM(1,2));
	set(S.dcm13,'string',num2str(DCM(1,3),'%8.6f'),'value',DCM(1,3));
	set(S.dcm21,'string',num2str(DCM(2,1),'%8.6f'),'value',DCM(2,1));
	set(S.dcm22,'string',num2str(DCM(2,2),'%8.6f'),'value',DCM(2,2));
	set(S.dcm23,'string',num2str(DCM(2,3),'%8.6f'),'value',DCM(2,3));
	set(S.dcm31,'string',num2str(DCM(3,1),'%8.6f'),'value',DCM(3,1));
	set(S.dcm32,'string',num2str(DCM(3,2),'%8.6f'),'value',DCM(3,2));
	set(S.dcm33,'string',num2str(DCM(3,3),'%8.6f'),'value',DCM(3,3));
	set(S.q1,'string',num2str(Q(1),'%8.6f'),'value',Q(1));
	set(S.q2,'string',num2str(Q(2),'%8.6f'),'value',Q(2));
	set(S.q3,'string',num2str(Q(3),'%8.6f'),'value',Q(3));
	set(S.q4,'string',num2str(Q(4),'%8.6f'),'value',Q(4));
	set(S.ep1,'string',num2str(EV(1),'%8.6f'),'value',EV(1));
	set(S.ep2,'string',num2str(EV(2),'%8.6f'),'value',EV(2));
	set(S.ep3,'string',num2str(EV(3),'%8.6f'),'value',EV(3));
	set(S.ep4,'string',num2str(EV(4),'%8.6f'),'value',EV(4));
	set(S.dialog,'string',errorstring);
    if ~strcmpi(errorstring(1,1:5),'error')
        S=plotter(S);
    end
end

function []=spincalcdcm(varargin)
    S=varargin{3};
	dcm11=get(S.dcm11,'value');
	dcm12=get(S.dcm12,'value');
	dcm13=get(S.dcm13,'value');
    dcm21=get(S.dcm21,'value');
	dcm22=get(S.dcm22,'value');
	dcm23=get(S.dcm23,'value');
    dcm31=get(S.dcm31,'value');
	dcm32=get(S.dcm32,'value');
	dcm33=get(S.dcm33,'value');
    dcm=[dcm11,dcm12,dcm13;dcm21,dcm22,dcm23;dcm31,dcm32,dcm33];
	EA_rotation_order_index=get(S.eamenu,'value');
	EA_rotations={'EA123';'EA132';'EA213';'EA231';'EA312';'EA321';'EA121';'EA131';'EA212';'EA232';'EA313';'EA323'};
	EA_rotation_order=EA_rotations{EA_rotation_order_index,:};
	[EA,errorstring]=spincalcmod(['DCMto',EA_rotation_order],dcm,1e-5,1);
	[EV,errorstring]=spincalcmod('DCMtoEV',dcm,1e-10,1);
	[Q,errorstring]=spincalcmod('DCMtoQ',dcm,1e-10,1);
    set(S.q1,'string',num2str(Q(1),'%8.6f'),'value',Q(1));
	set(S.q2,'string',num2str(Q(2),'%8.6f'),'value',Q(2));
	set(S.q3,'string',num2str(Q(3),'%8.6f'),'value',Q(3));
	set(S.q4,'string',num2str(Q(4),'%8.6f'),'value',Q(4));
	set(S.ep1,'string',num2str(EV(1),'%8.6f'),'value',EV(1));
	set(S.ep2,'string',num2str(EV(2),'%8.6f'),'value',EV(2));
	set(S.ep3,'string',num2str(EV(3),'%8.6f'),'value',EV(3));
	set(S.ep4,'string',num2str(EV(4),'%8.6f'),'value',EV(4));
    set(S.ea1,'string',num2str(EA(1),'%8.6f'),'value',EA(1));
	set(S.ea2,'string',num2str(EA(2),'%8.6f'),'value',EA(2));
	set(S.ea3,'string',num2str(EA(3),'%8.6f'),'value',EA(3));
	set(S.dialog,'string',errorstring);
    if ~strcmpi(errorstring(1,1:5),'error')
        S=plotter(S);
    end
end

function []=spincalcq(varargin)
	S=varargin{3};
	q1=get(S.q1,'value');
	q2=get(S.q2,'value');
	q3=get(S.q3,'value');
	q4=get(S.q4,'value');
	EA_rotation_order_index=get(S.eamenu,'value');
	EA_rotations={'EA123';'EA132';'EA213';'EA231';'EA312';'EA321';'EA121';'EA131';'EA212';'EA232';'EA313';'EA323'};
	EA_rotation_order=EA_rotations{EA_rotation_order_index,:};
	[DCM,errorstring]=spincalcmod('QtoDCM',[q1,q2,q3,q4],eps,1);	
	[EV,errorstring]=spincalcmod('QtoEV',[q1,q2,q3,q4],eps,1);
	[EA,errorstring]=spincalcmod(['Qto',EA_rotation_order],[q1,q2,q3,q4],eps,1);
	set(S.dcm11,'string',num2str(DCM(1,1),'%8.6f'),'value',DCM(1,1));
	set(S.dcm12,'string',num2str(DCM(1,2),'%8.6f'),'value',DCM(1,2));
	set(S.dcm13,'string',num2str(DCM(1,3),'%8.6f'),'value',DCM(1,3));
	set(S.dcm21,'string',num2str(DCM(2,1),'%8.6f'),'value',DCM(2,1));
	set(S.dcm22,'string',num2str(DCM(2,2),'%8.6f'),'value',DCM(2,2));
	set(S.dcm23,'string',num2str(DCM(2,3),'%8.6f'),'value',DCM(2,3));
	set(S.dcm31,'string',num2str(DCM(3,1),'%8.6f'),'value',DCM(3,1));
	set(S.dcm32,'string',num2str(DCM(3,2),'%8.6f'),'value',DCM(3,2));
	set(S.dcm33,'string',num2str(DCM(3,3),'%8.6f'),'value',DCM(3,3));
	set(S.ep1,'string',num2str(EV(1),'%8.6f'),'value',EV(1));
	set(S.ep2,'string',num2str(EV(2),'%8.6f'),'value',EV(2));
	set(S.ep3,'string',num2str(EV(3),'%8.6f'),'value',EV(3));
	set(S.ep4,'string',num2str(EV(4),'%8.6f'),'value',EV(4));
	set(S.ea1,'string',num2str(EA(1),'%8.6f'),'value',EA(1));
	set(S.ea2,'string',num2str(EA(2),'%8.6f'),'value',EA(2));
	set(S.ea3,'string',num2str(EA(3),'%8.6f'),'value',EA(3));
	set(S.dialog,'string',errorstring);
    if ~strcmpi(errorstring(1,1:5),'error')
        S=plotter(S);
    end
end

function []=spincalcep(varargin)
	S=varargin{3};
	m1=get(S.ep1,'value');
	m2=get(S.ep2,'value');
	m3=get(S.ep3,'value');
	mu=get(S.ep4,'value');
	EA_rotation_order_index=get(S.eamenu,'value');
	EA_rotations={'EA123';'EA132';'EA213';'EA231';'EA312';'EA321';'EA121';'EA131';'EA212';'EA232';'EA313';'EA323'};
	EA_rotation_order=EA_rotations{EA_rotation_order_index,:};
	[DCM,errorstring]=spincalcmod('EVtoDCM',[m1,m2,m3,mu],eps,1);	
	[Q,errorstring]=spincalcmod('EVtoQ',[m1,m2,m3,mu],eps,1);
	[EA,errorstring]=spincalcmod(['EVto',EA_rotation_order],[m1,m2,m3,mu],eps,1);
	set(S.dcm11,'string',num2str(DCM(1,1),'%8.6f'),'value',DCM(1,1));
	set(S.dcm12,'string',num2str(DCM(1,2),'%8.6f'),'value',DCM(1,2));
	set(S.dcm13,'string',num2str(DCM(1,3),'%8.6f'),'value',DCM(1,3));
	set(S.dcm21,'string',num2str(DCM(2,1),'%8.6f'),'value',DCM(2,1));
	set(S.dcm22,'string',num2str(DCM(2,2),'%8.6f'),'value',DCM(2,2));
	set(S.dcm23,'string',num2str(DCM(2,3),'%8.6f'),'value',DCM(2,3));
	set(S.dcm31,'string',num2str(DCM(3,1),'%8.6f'),'value',DCM(3,1));
	set(S.dcm32,'string',num2str(DCM(3,2),'%8.6f'),'value',DCM(3,2));
	set(S.dcm33,'string',num2str(DCM(3,3),'%8.6f'),'value',DCM(3,3));
	set(S.q1,'string',num2str(Q(1),'%8.6f'),'value',Q(1));
	set(S.q2,'string',num2str(Q(2),'%8.6f'),'value',Q(2));
	set(S.q3,'string',num2str(Q(3),'%8.6f'),'value',Q(3));
	set(S.q4,'string',num2str(Q(4),'%8.6f'),'value',Q(4));
	set(S.ea1,'string',num2str(EA(1),'%8.6f'),'value',EA(1));
	set(S.ea2,'string',num2str(EA(2),'%8.6f'),'value',EA(2));
	set(S.ea3,'string',num2str(EA(3),'%8.6f'),'value',EA(3));
	set(S.dialog,'string',errorstring);
    if ~strcmpi(errorstring(1,1:5),'error')
        S=plotter(S);
    end
end

function [OUTPUT,errorstring]=spincalcmod(CONVERSION,INPUT,tol,ichk)
%Function for the conversion of one rotation input type to desired output.
%Supported conversion input/output types are as follows:
%   1: Q        Rotation Quaternions
%   2: EV       Euler Vector and rotation angle (degrees)
%   3: DCM      Orthogonal DCM Rotation Matrix
%   4: EA###    Euler angles (12 possible sets) (degrees)
%
%Author: John Fuller
%National Institute of Aerospace
%Hampton, VA 23666
%John.Fuller@nianet.org
%
%Version 1.3
%June 30th, 2009
%
%Version 1.3 updates
%   SpinCalc now detects when input data is too close to Euler singularity, if user is choosing
%   Euler angle output. Prohibits output if middle angle is within 0.1 degree of singularity value.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                OUTPUT=spincalcmod(CONVERSION,INPUT,tol,ichk)
%Inputs:
%CONVERSION - Single string value that dictates the type of desired
%             conversion.  The conversion strings are listed below.
%
%   'DCMtoEA###'  'DCMtoEV'    'DCMtoQ'       **for cases that involve   
%   'EA###toDCM'  'EA###toEV'  'EA###toQ'       euler angles, ### should be
%   'EVtoDCM'     'EVtoEA###'  'EVtoQ'          replaced with the proper 
%   'QtoDCM'      'QtoEA###'   'QtoEV'          order desired.  EA321 would
%   'EA###toEA###'                              be Z(yaw)-Y(pitch)-X(roll).
%
%INPUT - matrix or vector that corresponds to the first entry in the
%        CONVERSION string, formatted as follows:
%
%        DCM - 3x3xN multidimensional matrix which pre-multiplies by a coordinate
%              frame vector to rotate it to the desired new frame.
%
%        EA### - [psi,theta,phi] (Nx3) row vector list dictating to the first angle
%                rotation (psi), the second (theta), and third (phi) (DEGREES)
%
%        EV - [m1,m2,m3,MU] (Nx4) row vector list dictating the components of euler
%             rotation vector (original coordinate frame) and the Euler 
%             rotation angle about that vector (MU) (DEGREES)
%
%        Q - [q1,q2,q3,q4] (Nx4) row vector list defining quaternion of
%            rotation.  q4 = cos(MU/2) where MU is Euler rotation angle
%
%tol - tolerance value
%ichk - 0 disables warning flags
%          1 enables warning flags (near singularities)
%**NOTE: N corresponds to multiple orientations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Output:
%OUTPUT - matrix or vector corresponding to the second entry in the
%         CONVERSION input string, formatted as shown above.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Pre-processer to determine type of conversion from CONVERSION string input
%Types are numbered as follows:
%Q=1   EV=2   DCM=3   EA=4
i_type=strfind(lower(CONVERSION),'to');
length=size(CONVERSION,2);
error_flag=0;
errorstring='No errors found.';
if length>12 || length<4,   %no CONVERSION string can be shorter than 4 or longer than 12 chars
    error('Error: Invalid entry for CONVERSION input string');
end
o_type=length-i_type;
if i_type<5,
    i_type=i_type-1;
else
    i_type=i_type-2;
end
if o_type<5,
    o_type=o_type-1;
else
    o_type=o_type-2;
end
TYPES=cell(1,4);
TYPES{1,1}='Q'; TYPES{1,2}='EV'; TYPES{1,3}='DCM'; TYPES{1,4}='EA';
INPUT_TYPE=TYPES{1,i_type};
OUTPUT_TYPE=TYPES{1,o_type};
clear TYPES
%Confirm input as compared to program interpretation
if i_type~=4 && o_type~=4,  %if input/output are NOT Euler angles
    CC=[INPUT_TYPE,'to',OUTPUT_TYPE];
    if strcmpi(CONVERSION,CC)==0;
        error('Error: Invalid entry for CONVERSION input string');        
    end
else
    if i_type==4,   %if input type is Euler angles, determine the order of rotations
        EULER_order_in=str2double(CONVERSION(1,3:5));
        rot_1_in=floor(EULER_order_in/100);     %first rotation
        rot_2_in=floor((EULER_order_in-rot_1_in*100)/10);   %second rotation
        rot_3_in=(EULER_order_in-rot_1_in*100-rot_2_in*10);   %third rotation
        if rot_1_in<1 || rot_2_in<1 || rot_3_in<1 || rot_1_in>3 || rot_2_in>3 || rot_3_in>3,
            error('Error: Invalid input Euler angle order type (conversion string).');  %check that all orders are between 1 and 3            
        elseif rot_1_in==rot_2_in || rot_2_in==rot_3_in,
            error('Error: Invalid input Euler angle order type (conversion string).');  %check that no 2 consecutive orders are equal (invalid)           
        end
        %check input dimensions to be 1x3x1
        if size(INPUT,2)~=3 || size(INPUT,3)~=1
            error('Error: Input euler angle data vector is not Nx3')            
        end
        %identify singularities        
        if rot_1_in==rot_3_in, %Type 2 rotation (first and third rotations about same axis)          
            if INPUT(:,2)<=zeros(size(INPUT,1),1) | INPUT(:,2)>=180*ones(size(INPUT,1),1),  %#ok<OR2> %confirm second angle within range
                %error('Error: Second input Euler angle(s) outside 0 to 180 degree range')
				errorstring='Error: Second input Euler angle(s) outside 0 to 180 degree range';
				error_flag=1;
            elseif abs(INPUT(:,2))<2*ones(size(INPUT,1),1) | abs(INPUT(:,2))>178*ones(size(INPUT,1),1),  %#ok<OR2> %check for singularity
                if ichk==1,
                    %errordlg('Warning: Input Euler angle rotation(s) near a singularity.               Second angle near 0 or 180 degrees.')
					errorstring={'Warning: Input Euler angle rotation(s) near a singularity.';'Second angle near 0 or 180 degrees.'};
                end
            end
        else    %Type 1 rotation (all rotations about each of three axes)
            if abs(INPUT(:,2))>=90*ones(size(INPUT,1),1), %confirm second angle within range
                %error('Error: Second input Euler angle(s) outside -90 to 90 degree range')
				errorstring='Error: Second input Euler angle(s) outside -90 to 90 degree range';
				error_flag=1;
            elseif abs(INPUT(:,2))>88*ones(size(INPUT,1),1),  %check for singularity
                if ichk==1, %#ok<ALIGN>
                    %errordlg('Warning: Input Euler angle(s) rotation near a singularity.               Second angle near -90 or 90 degrees.')
					errorstring={'Warning: Input Euler angle(s) rotation near a singularity.';'Second angle near -90 or 90 degrees.'};
				end
            end
        end    
    end
    if o_type==4,   %if output type is Euler angles, determine order of rotations
        EULER_order_out=str2double(CONVERSION(1,length-2:length));
        rot_1_out=floor(EULER_order_out/100);   %first rotation
        rot_2_out=floor((EULER_order_out-rot_1_out*100)/10);    %second rotation
        rot_3_out=(EULER_order_out-rot_1_out*100-rot_2_out*10); %third rotation
        if rot_1_out<1 || rot_2_out<1 || rot_3_out<1 || rot_1_out>3 || rot_2_out>3 || rot_3_out>3,
            error('Error: Invalid output Euler angle order type (conversion string).'); %check that all orders are between 1 and 3           
        elseif rot_1_out==rot_2_out || rot_2_out==rot_3_out,
            error('Error: Invalid output Euler angle order type (conversion string).'); %check that no 2 consecutive orders are equal
        end
    end        
    if i_type==4 && o_type~=4,  %if input are euler angles but not output
        CC=['EA',num2str(EULER_order_in),'to',OUTPUT_TYPE]; %construct program conversion string for checking against user input
    elseif o_type==4 && i_type~=4,  %if output are euler angles but not input
        CC=[INPUT_TYPE,'to','EA',num2str(EULER_order_out)]; %construct program conversion string for checking against user input
    elseif i_type==4 && o_type==4,  %if both input and output are euler angles
        CC=['EA',num2str(EULER_order_in),'to','EA',num2str(EULER_order_out)];   %construct program conversion string
    end
    if strcmpi(CONVERSION,CC)==0; %check program conversion string against user input to confirm the conversion command
        error('Error: Invalid entry for CONVERSION input string');
    end
end
clear i_type o_type CC

%From the input, determine the quaternions that uniquely describe the
%rotation prescribed by that input.  The output will be calculated in the
%second portion of the code from these quaternions.
switch INPUT_TYPE
    case 'DCM'
        if size(INPUT,1)~=3 || size(INPUT,2)~=3  %check DCM dimensions
            error('Error: DCM matrix is not 3x3xN');           
        end
        N=size(INPUT,3);    %number of orientations
        %Check if matrix is indeed orthogonal
        perturbed=NaN(3,3,N);
        DCM_flag=0;
        for ii=1:N,
            perturbed(:,:,ii)=abs(INPUT(:,:,ii)*INPUT(:,:,ii)'-eye(3)); %perturbed array shows difference between DCM*DCM' and I
            if abs(det(INPUT(:,:,ii))-1)>tol, %if determinant is off by one more than tol, user is warned.
                if ichk==1,
                    DCM_flag=1;
                end
            end
            if abs(det(INPUT(:,:,ii))+1)<0.05, %if determinant is near -1, DCM is improper
                %error('Error: Input DCM(s) improper');
				errorstring='Error: Input DCM(s) improper.';
				error_flag=1;
				break
            end
            if DCM_flag==1,
                %errordlg('Warning: Input DCM matrix determinant(s) off from 1 by more than tolerance.')
				errorstring='Warning: Input DCM matrix determinant(s) off from 1 by more than tolerance.';
            end
        end
        DCM_flag=0;
        if ichk==1,
            for kk=1:N,
                for ii=1:3,
                    for jj=1:3,
                        if perturbed(ii,jj,kk)>tol,   %if any difference is larger than tol, user is warned.
                            DCM_flag=1;
                        end
                    end
                end
            end
            if DCM_flag==1,
                %fprintf('Warning: Input DCM(s) matrix not orthogonal to precision tolerance.')
				errorstring='Warning: Input DCM(s) matrix not orthogonal to precision tolerance.';
            end
        end       
        clear perturbed DCM_flag   
        Q=NaN(4,N);
        for ii=1:N,
            denom=NaN(4,1);
            denom(1)=0.5*sqrt(1+INPUT(1,1,ii)-INPUT(2,2,ii)-INPUT(3,3,ii));
            denom(2)=0.5*sqrt(1-INPUT(1,1,ii)+INPUT(2,2,ii)-INPUT(3,3,ii));
            denom(3)=0.5*sqrt(1-INPUT(1,1,ii)-INPUT(2,2,ii)+INPUT(3,3,ii));
            denom(4)=0.5*sqrt(1+INPUT(1,1,ii)+INPUT(2,2,ii)+INPUT(3,3,ii));        
            %determine which Q equations maximize denominator
            switch find(denom==max(denom),1,'first')  %determines max value of qtests to put in denominator
                case 1
                    Q(1,ii)=denom(1);
                    Q(2,ii)=(INPUT(1,2,ii)+INPUT(2,1,ii))/(4*Q(1,ii));
                    Q(3,ii)=(INPUT(1,3,ii)+INPUT(3,1,ii))/(4*Q(1,ii));
                    Q(4,ii)=(INPUT(2,3,ii)-INPUT(3,2,ii))/(4*Q(1,ii));
                case 2
                    Q(2,ii)=denom(2);
                    Q(1,ii)=(INPUT(1,2,ii)+INPUT(2,1,ii))/(4*Q(2,ii));
                    Q(3,ii)=(INPUT(2,3,ii)+INPUT(3,2,ii))/(4*Q(2,ii));
                    Q(4,ii)=(INPUT(3,1,ii)-INPUT(1,3,ii))/(4*Q(2,ii));
                case 3
                    Q(3,ii)=denom(3);
                    Q(1,ii)=(INPUT(1,3,ii)+INPUT(3,1,ii))/(4*Q(3,ii));
                    Q(2,ii)=(INPUT(2,3,ii)+INPUT(3,2,ii))/(4*Q(3,ii));
                    Q(4,ii)=(INPUT(1,2,ii)-INPUT(2,1,ii))/(4*Q(3,ii));
                case 4
                    Q(4,ii)=denom(4);
                    Q(1,ii)=(INPUT(2,3,ii)-INPUT(3,2,ii))/(4*Q(4,ii));
                    Q(2,ii)=(INPUT(3,1,ii)-INPUT(1,3,ii))/(4*Q(4,ii));
                    Q(3,ii)=(INPUT(1,2,ii)-INPUT(2,1,ii))/(4*Q(4,ii));
            end
        end
        Q=Q';
        clear denom
    case 'EV'  %Euler Vector Input Type
        if size(INPUT,2)~=4 || size(INPUT,3)~=1   %check dimensions
            error('Error: Input euler vector and rotation data matrix is not Nx4')            
        end
        N=size(INPUT,1);
        MU=INPUT(:,4)*pi/180;  %assign mu name for clarity
        if abs(sqrt(INPUT(:,1).^2+INPUT(:,2).^2+INPUT(:,3).^2)-ones(N,1))>tol*ones(N,1),  %check that input m's constitute unit vector
            %error('Input euler vector(s) components do not constitute a unit vector')
			errorstring='Error: Input euler vector(s) components do not constitute a unit vector.';
			error_flag=1;
        end
        if MU<-2*pi*ones(N,1) || MU>2*pi*ones(N,1), %check if rotation about euler vector is between 0 and 360
            %error('Input euler rotation angle(s) not between -360 and 360 degrees')
			errorstring='Error: Input mu rotation angle(s) not between -360 and 360 degrees.';
			error_flag=1;
        end
        Q=[INPUT(:,1).*sin(MU/2),INPUT(:,2).*sin(MU/2),INPUT(:,3).*sin(MU/2),cos(MU/2)];   %quaternion
        clear m1 m2 m3 MU
    case 'EA'        
        psi=INPUT(:,1)*pi/180;  theta=INPUT(:,2)*pi/180;  phi=INPUT(:,3)*pi/180;
        N=size(INPUT,1);    %number of orientations
        %Pre-calculate cosines and sines of the half-angles for conversion.
        c1=cos(psi./2); c2=cos(theta./2); c3=cos(phi./2);
        s1=sin(psi./2); s2=sin(theta./2); s3=sin(phi./2);
        c13=cos((psi+phi)./2);  s13=sin((psi+phi)./2);
        c1_3=cos((psi-phi)./2);  s1_3=sin((psi-phi)./2);
        c3_1=cos((phi-psi)./2);  s3_1=sin((phi-psi)./2);
        if EULER_order_in==121,
            Q=[c2.*s13,s2.*c1_3,s2.*s1_3,c2.*c13];
        elseif EULER_order_in==232,
            Q=[s2.*s1_3,c2.*s13,s2.*c1_3,c2.*c13];
        elseif EULER_order_in==313;
            Q=[s2.*c1_3,s2.*s1_3,c2.*s13,c2.*c13];
        elseif EULER_order_in==131,
            Q=[c2.*s13,s2.*s3_1,s2.*c3_1,c2.*c13];
        elseif EULER_order_in==212,
            Q=[s2.*c3_1,c2.*s13,s2.*s3_1,c2.*c13];
        elseif EULER_order_in==323,
            Q=[s2.*s3_1,s2.*c3_1,c2.*s13,c2.*c13];
        elseif EULER_order_in==123,
            Q=[s1.*c2.*c3+c1.*s2.*s3,c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*s3+s1.*s2.*c3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==231,
            Q=[c1.*c2.*s3+s1.*s2.*c3,s1.*c2.*c3+c1.*s2.*s3,c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==312,
            Q=[c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*s3+s1.*s2.*c3,s1.*c2.*c3+c1.*s2.*s3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==132,
            Q=[s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*s3-s1.*s2.*c3,c1.*s2.*c3+s1.*c2.*s3,c1.*c2.*c3+s1.*s2.*s3];
        elseif EULER_order_in==213,
            Q=[c1.*s2.*c3+s1.*c2.*s3,s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*s3-s1.*s2.*c3,c1.*c2.*c3+s1.*s2.*s3];
        elseif EULER_order_in==321,
            Q=[c1.*c2.*s3-s1.*s2.*c3,c1.*s2.*c3+s1.*c2.*s3,s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*c3+s1.*s2.*s3];
        else
            error('Error: Invalid input Euler angle order type (conversion string)');            
        end
        clear c1 s1 c2 s2 c3 s3 c13 s13 c1_3 s1_3 c3_1 s3_1 psi theta phi
    case 'Q'
        if size(INPUT,2)~=4 || size(INPUT,3)~=1
            error('Error: Input quaternion matrix is not Nx4');            
        end
        N=size(INPUT,1);    %number of orientations 
        if ichk==1,
            if abs(sqrt(INPUT(:,1).^2+INPUT(:,2).^2+INPUT(:,3).^2+INPUT(:,4).^2)-ones(N,1))>tol*ones(N,1)
                %errordlg('Warning: Input quaternion norm(s) deviate(s) from unity by more than tolerance')
				errorstring='Warning: Input quaternion norm(s) deviate(s) from unity by more than tolerance';
            end 
        end
        Q=INPUT;
end
clear INPUT INPUT_TYPE EULER_order_in

%Normalize quaternions in case of deviation from unity.  User has already
%been warned of deviation.
Qnorms=sqrt(sum(Q.*Q,2));
Q=[Q(:,1)./Qnorms,Q(:,2)./Qnorms,Q(:,3)./Qnorms,Q(:,4)./Qnorms];

switch OUTPUT_TYPE
	case 'DCM'
		OUTPUT=NaN(3,3);
	case 'EV'
		OUTPUT=NaN(1,4);
	case 'Q'
		OUTPUT=NaN(1,4);
	case 'EA'
		OUTPUT=NaN(1,3);
end

if error_flag==0
	switch OUTPUT_TYPE
		case 'DCM'
			Q=reshape(Q',1,4,N);
			OUTPUT=[Q(1,1,:).^2-Q(1,2,:).^2-Q(1,3,:).^2+Q(1,4,:).^2,2*(Q(1,1,:).*Q(1,2,:)+Q(1,3,:).*Q(1,4,:)),2*(Q(1,1,:).*Q(1,3,:)-Q(1,2,:).*Q(1,4,:));
					2*(Q(1,1,:).*Q(1,2,:)-Q(1,3,:).*Q(1,4,:)),-Q(1,1,:).^2+Q(1,2,:).^2-Q(1,3,:).^2+Q(1,4,:).^2,2*(Q(1,2,:).*Q(1,3,:)+Q(1,1,:).*Q(1,4,:));
					2*(Q(1,1,:).*Q(1,3,:)+Q(1,2,:).*Q(1,4,:)),2*(Q(1,2,:).*Q(1,3,:)-Q(1,1,:).*Q(1,4,:)),-Q(1,1,:).^2-Q(1,2,:).^2+Q(1,3,:).^2+Q(1,4,:).^2];
		case 'EV'
			MU=2*atan2(sqrt(sum(Q(:,1:3).*Q(:,1:3),2)),Q(:,4));
			if sin(MU/2)~=zeros(N,1),
				OUTPUT=[Q(:,1)./sin(MU/2),Q(:,2)./sin(MU/2),Q(:,3)./sin(MU/2),MU*180/pi];
			else
				OUTPUT=NaN(N,4);
				for ii=1:N,
					if sin(MU(ii,1)/2)~=0,
						OUTPUT(ii,1:4)=[Q(ii,1)/sin(MU(ii,1)/2),Q(ii,2)/sin(MU(ii,1)/2),Q(ii,3)/sin(MU(ii,1)/2),MU(ii,1)*180/pi];
					else
						OUTPUT(ii,1:4)=[1,0,0,MU(ii,1)*180/pi];
					end
				end
			end
		case 'Q'
			OUTPUT=Q;
		case 'EA'
			if EULER_order_out==121,
				psi=atan2((Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
				theta=acos(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2);
				phi=atan2((Q(:,1).*Q(:,2)-Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
		  Euler_type=2;
			elseif EULER_order_out==232;
				psi=atan2((Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
				theta=acos(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2);
				phi=atan2((Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)),(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
		  Euler_type=2;
			elseif EULER_order_out==313;
				psi=atan2((Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
				theta=acos(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2);
				phi=atan2((Q(:,1).*Q(:,3)-Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
		  Euler_type=2;
			elseif EULER_order_out==131;
				psi=atan2((Q(:,1).*Q(:,3)-Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
				theta=acos(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2);
				phi=atan2((Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
		  Euler_type=2;
			elseif EULER_order_out==212;
				psi=atan2((Q(:,1).*Q(:,2)-Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
				theta=acos(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2);
				phi=atan2((Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
		  Euler_type=2;
			elseif EULER_order_out==323;
				psi=atan2((Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)),(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
				theta=acos(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2);
				phi=atan2((Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
		  Euler_type=2;
			elseif EULER_order_out==123;
				psi=atan2(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
				theta=asin(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
				phi=atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
		  Euler_type=1;
			elseif EULER_order_out==231;
				psi=atan2(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
				theta=asin(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
				phi=atan2(2.*(Q(:,1).*Q(:,4)-Q(:,3).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
		  Euler_type=1;
			elseif EULER_order_out==312;
				psi=atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
				theta=asin(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
				phi=atan2(2.*(Q(:,2).*Q(:,4)-Q(:,3).*Q(:,1)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
		  Euler_type=1;
			elseif EULER_order_out==132;
				psi=atan2(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
				theta=asin(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
				phi=atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
		  Euler_type=1;
			elseif EULER_order_out==213;
				psi=atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
				theta=asin(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
				phi=atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
		  Euler_type=1;
			elseif EULER_order_out==321;
				psi=atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
				theta=asin(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
				phi=atan2(2.*(Q(:,1).*Q(:,4)+Q(:,3).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
		  Euler_type=1;
			else
				error('Error: Invalid output Euler angle order type (conversion string).');           
			end
			if(isreal([psi,theta,phi]))==0,
				%error('Error: Unreal Euler output.  Input resides too close to singularity.  Please choose different output type.')
				errorstring={'Error: Unreal Euler Angle output.  Input resides too close to singularity.';'Please choose different output type.'};
				OUTPUT=[NaN,NaN,NaN];
			end
			OUTPUT=mod([psi,theta,phi]*180/pi,360);  %deg
			if Euler_type==1, %#ok<ALIGN>
				sing_chk=find(abs(theta)*180/pi>89.9);
				sing_chk=sort(sing_chk(sing_chk>0));
				if size(sing_chk,1)>=1,
					%error('Error: Input rotation #%s resides too close to Type 1 Euler singularity.\nType 1 Euler singularity occurs when second angle is -90 or 90 degrees.\nPlease choose different output type.',num2str(sing_chk(1,1)));
					errorstring={'Error: Input rotation resides too close to Type 1 Euler singularity.';'Type 1 Euler singularity occurs when second angle is -90 or 90 degrees.';'Please choose a different EA order.'};
					OUTPUT=[NaN,NaN,NaN];
				end
			elseif Euler_type==2,
				sing_chk=[find(abs(theta*180/pi)<0.1);find(abs(theta*180/pi-180)<0.1);find(abs(theta*180/pi-360))<0.1];
				sing_chk=sort(sing_chk(sing_chk>0));
				if size(sing_chk,1)>=1, %#ok<ALIGN>
					%error('Error: Input rotation #%s resides too close to Type 2 Euler singularity.\nType 2 Euler singularity occurs when second angle is 0 or 180 degrees.\nPlease choose different output type.',num2str(sing_chk(1,1)));
					errorstring={'Error: Input rotation resides too close to Type 2 Euler singularity.';'Type 2 Euler singularity occurs when second angle is 0 or 180 degrees.';'Please choose a different EA order.'};
					OUTPUT=[NaN,NaN,NaN];
                end
            end
            %Modified output Euler angles to be between -180 and 180
            temp=OUTPUT(:,1);
            temp(temp>180)=temp(temp>180)-360;
            OUTPUT(:,1)=temp;
            temp=OUTPUT(:,2);
            temp(temp>180)=temp(temp>180)-360;
            OUTPUT(:,2)=temp;
            temp=OUTPUT(:,3);
            temp(temp>180)=temp(temp>180)-360;
            OUTPUT(:,3)=temp;
	end
end
OUTPUT(abs(OUTPUT)<1e-14)=0;
end

function []=closegui(varargin)
    S=varargin{3};
    close(S.fh)
end

function []=assignep1(varargin)
    S=varargin{3};
    set(S.ep1,'value',str2double(get(S.ep1,'string')));
end

function []=assignep2(varargin)
    S=varargin{3};
    set(S.ep2,'value',str2double(get(S.ep2,'string')));
end

function []=assignep3(varargin)
    S=varargin{3};
    set(S.ep3,'value',str2double(get(S.ep3,'string')));
end

function []=assignea1(varargin)
    S=varargin{3};
    set(S.ea1,'value',str2double(get(S.ea1,'string')));
end

function []=assignea2(varargin)
    S=varargin{3};
    set(S.ea2,'value',str2double(get(S.ea2,'string')));
end

function []=assignea3(varargin)
    S=varargin{3};
    set(S.ea3,'value',str2double(get(S.ea3,'string')));
end

function []=assignep4(varargin)
    S=varargin{3};
    set(S.ep4,'value',str2double(get(S.ep4,'string')));
end

function []=assignq1(varargin)
    S=varargin{3};
    set(S.q1,'value',str2double(get(S.q1,'string')));
end

function []=assignq2(varargin)
    S=varargin{3};
    set(S.q2,'value',str2double(get(S.q2,'string')));
end

function []=assignq3(varargin)
    S=varargin{3};
    set(S.q3,'value',str2double(get(S.q3,'string')));
end

function []=assignq4(varargin)
    S=varargin{3};
    set(S.q4,'value',str2double(get(S.q4,'string')));
end

function []=assigndcm11(varargin)
    S=varargin{3};
    set(S.dcm11,'value',str2double(get(S.dcm11,'string')));
end

function []=assigndcm12(varargin)
    S=varargin{3};
    set(S.dcm12,'value',str2double(get(S.dcm12,'string')));
end

function []=assigndcm13(varargin)
    S=varargin{3};
    set(S.dcm13,'value',str2double(get(S.dcm13,'string')));
end

function []=assigndcm21(varargin)
    S=varargin{3};
    set(S.dcm21,'value',str2double(get(S.dcm21,'string')));
end

function []=assigndcm22(varargin)
    S=varargin{3};
    set(S.dcm22,'value',str2double(get(S.dcm22,'string')));
end

function []=assigndcm23(varargin)
    S=varargin{3};
    set(S.dcm23,'value',str2double(get(S.dcm23,'string')));
end

function []=assigndcm31(varargin)
    S=varargin{3};
    set(S.dcm31,'value',str2double(get(S.dcm31,'string')));
end

function []=assigndcm32(varargin)
    S=varargin{3};
    set(S.dcm32,'value',str2double(get(S.dcm32,'string')));
end

function []=assigndcm33(varargin)
    S=varargin{3};
    set(S.dcm33,'value',str2double(get(S.dcm33,'string')));
end

function [hout,ax_out] = uibutton(varargin)
%uibutton: Create pushbutton with more flexible labeling than uicontrol.
% Usage:
%   uibutton accepts all the same arguments as uicontrol except for the
%   following property changes:
%
%     Property      Values
%     -----------   ------------------------------------------------------
%     Style         'pushbutton', 'togglebutton' or 'text', default =
%                   'pushbutton'.
%     String        Same as for text() including cell array of strings and
%                   TeX or LaTeX interpretation.
%     Interpreter   'tex', 'latex' or 'none', default = default for text()
%     Rotation      text rotation angle, default = 0
%
% Syntax:
%   handle = uibutton('PropertyName',PropertyValue,...)
%   handle = uibutton(parent,'PropertyName',PropertyValue,...)
%   [text_obj,axes_handle] = uibutton('Style','text',...
%       'PropertyName',PropertyValue,...)
%
% uibutton creates a temporary axes and text object containing the text to
% be displayed, captures the axes as an image, deletes the axes and then
% displays the image on the uicontrol.  The handle to the uicontrol is
% returned.  If you pass in a handle to an existing uicontol as the first
% argument then uibutton will use that uicontrol and not create a new one.
%
% If the Style is set to 'text' then the axes object is not deleted and the
% text object handle is returned (as well as the handle to the axes in a
% second output argument).
%
% See also UICONTROL.

% Version: 1.8, 10 March 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Detect if first argument is a uicontrol handle.
keep_handle = false;
if nargin > 0
	h = varargin{1};
	if isscalar(h) && ishandle(h) && strcmp(get(h,'Type'),'uicontrol')
		keep_handle = true;
		varargin(1) = [];
	end
end

% Parse arguments looking for 'Interpreter' property.  If found, note its
% value and then remove it from where it was found.
interp_value = get(0,'DefaultTextInterpreter');
rotation_value = get(0,'DefaultTextRotation');
arg = 1;
remove = [];
while arg <= length(varargin)
	v = varargin{arg};
	if isstruct(v)
		fn = fieldnames(v);
		for i = 1:length(fn)
			if strncmpi(fn{i},'interpreter',length(fn{i}))
				interp_value = v.(fn{i});
				v = rmfield(v,fn{i});
			elseif strncmpi(fn{i},'rotation',length(fn{i}))
				rotation_value = v.(fn{i});
				v = rmfield(v,fn{i});
			end
		end
		varargin{arg} = v;
		arg = arg + 1;
	elseif ischar(v)
		if strncmpi(v,'interpreter',length(v))
			interp_value = varargin{arg+1};
			remove = [remove,arg,arg+1]; %#ok<AGROW>
		elseif strncmpi(v,'rotation',length(v))
			rotation_value = varargin{arg+1};
			remove = [remove,arg,arg+1]; %#ok<AGROW>
		end
		arg = arg + 2;
	elseif arg == 1 && isscalar(v) && ishandle(v) && ...
			any(strcmp(get(h,'Type'),{'figure','uipanel'}))
		arg = arg + 1;
	else
		error('Invalid property or uicontrol parent.')
	end
end
varargin(remove) = [];

% Create uicontrol, get its properties then hide it.
if keep_handle
	set(h,varargin{:})
else
	h = uicontrol(varargin{:});
end
s = get(h);
if ~any(strcmp(s.Style,{'pushbutton','togglebutton','text'}))
	delete(h)
	error('''Style'' must be pushbutton, togglebutton or text.')
end
set(h,'Visible','off')

% Create axes.
parent = get(h,'Parent');
ax = axes('Parent',parent,...
	'Units',s.Units,...
	'Position',s.Position,...
	'XTick',[],'YTick',[],...
	'XColor',s.BackgroundColor,...
	'YColor',s.BackgroundColor,...
	'Box','on',...
	'Color',s.BackgroundColor);
% Adjust size of axes for best appearance.
set(ax,'Units','pixels')
pos = round(get(ax,'Position'));
if strcmp(s.Style,'text')
	set(ax,'Position',pos + [0 1 -1 -1])
else
	set(ax,'Position',pos + [4 4 -8 -8])
end
switch s.HorizontalAlignment
	case 'left'
		x = 0.0;
	case 'center'
		x = 0.5;
	case 'right'
		x = 1;
end
% Create text object.
text_obj = text('Parent',ax,...
	'Position',[x,0.5],...
	'String',s.String,...
	'Interpreter',interp_value,...
	'Rotation',rotation_value,...
	'HorizontalAlignment',s.HorizontalAlignment,...
	'VerticalAlignment','middle',...
	'FontName',s.FontName,...
	'FontSize',s.FontSize,...
	'FontAngle',s.FontAngle,...
	'FontWeight',s.FontWeight,...
	'Color',s.ForegroundColor);

% If we are creating something that looks like a text uicontrol then we're
% all done and we return the text object and axes handles rather than a
% uicontrol handle.
if strcmp(s.Style,'text')
	delete(h)
	if nargout
		hout = text_obj;
		ax_out = ax;
	end
	return
end

% Capture image of axes and then delete the axes.
frame = getframe(ax);
delete(ax)

% Build RGB image, set background pixels to NaN and put it in 'CData' for
% the uicontrol.
if isempty(frame.colormap)
	rgb = frame.cdata;
else
	rgb = reshape(frame.colormap(frame.cdata,:),[pos([4,3]),3]);
end
size_rgb = size(rgb);
rgb = double(rgb)/255;
back = repmat(permute(s.BackgroundColor,[1 3 2]),size_rgb(1:2));
isback = all(rgb == back,3);
rgb(repmat(isback,[1 1 3])) = NaN;
set(h,'CData',rgb,'String','','Visible',s.Visible)

% Assign output argument if necessary.
if nargout
	hout = h;
end
end