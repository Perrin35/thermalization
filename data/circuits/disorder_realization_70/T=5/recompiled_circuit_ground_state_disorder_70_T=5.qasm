OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.313247591257095) q[0];
sx q[0];
rz(3.91754099925096) q[0];
sx q[0];
rz(10.5269634485166) q[0];
rz(1.6232396364212) q[1];
sx q[1];
rz(2.01453879674012) q[1];
sx q[1];
rz(9.19934993087455) q[1];
cx q[1],q[0];
rz(0.0910873934626579) q[0];
sx q[0];
rz(3.84141811926896) q[0];
sx q[0];
rz(10.8989449500959) q[0];
rz(-1.36389768123627) q[2];
sx q[2];
rz(3.7676017006212) q[2];
sx q[2];
rz(10.7204563379209) q[2];
cx q[2],q[1];
rz(-0.0825833827257156) q[1];
sx q[1];
rz(4.34828522999818) q[1];
sx q[1];
rz(10.7081798076551) q[1];
rz(2.06462502479553) q[3];
sx q[3];
rz(4.40282371838624) q[3];
sx q[3];
rz(10.4260110616605) q[3];
cx q[3],q[2];
rz(0.483083575963974) q[2];
sx q[2];
rz(3.81566390593583) q[2];
sx q[2];
rz(9.40413675307437) q[2];
rz(0.189272999763489) q[3];
sx q[3];
rz(3.49116644461686) q[3];
sx q[3];
rz(11.0921983480374) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.780576884746552) q[0];
sx q[0];
rz(3.92266211112077) q[0];
sx q[0];
rz(10.4038549423139) q[0];
rz(-0.111152209341526) q[1];
sx q[1];
rz(3.87927553256089) q[1];
sx q[1];
rz(8.3637916803281) q[1];
cx q[1],q[0];
rz(1.84213352203369) q[0];
sx q[0];
rz(3.85028305848176) q[0];
sx q[0];
rz(10.243283545963) q[0];
rz(1.57681441307068) q[2];
sx q[2];
rz(3.00413470168645) q[2];
sx q[2];
rz(9.89203388094112) q[2];
cx q[2],q[1];
rz(2.09095215797424) q[1];
sx q[1];
rz(2.05558577378327) q[1];
sx q[1];
rz(9.65837713181182) q[1];
rz(1.95986974239349) q[3];
sx q[3];
rz(4.57923868496949) q[3];
sx q[3];
rz(7.87177107333347) q[3];
cx q[3],q[2];
rz(1.70791184902191) q[2];
sx q[2];
rz(4.96249893506105) q[2];
sx q[2];
rz(9.02672333120509) q[2];
rz(-1.71786439418793) q[3];
sx q[3];
rz(2.85456961591775) q[3];
sx q[3];
rz(9.9106013238351) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.61755359172821) q[0];
sx q[0];
rz(3.65312138398225) q[0];
sx q[0];
rz(10.4771133422773) q[0];
rz(-0.443183362483978) q[1];
sx q[1];
rz(4.10301503737504) q[1];
sx q[1];
rz(11.9124212026517) q[1];
cx q[1],q[0];
rz(0.674146473407745) q[0];
sx q[0];
rz(3.9735439141565) q[0];
sx q[0];
rz(9.25287426113292) q[0];
rz(0.623183608055115) q[2];
sx q[2];
rz(4.653548630076) q[2];
sx q[2];
rz(10.1544424056928) q[2];
cx q[2],q[1];
rz(-0.200507700443268) q[1];
sx q[1];
rz(4.294910581904) q[1];
sx q[1];
rz(11.7862293481748) q[1];
rz(-0.177860513329506) q[3];
sx q[3];
rz(4.486889513331) q[3];
sx q[3];
rz(10.0653862714688) q[3];
cx q[3],q[2];
rz(0.370864808559418) q[2];
sx q[2];
rz(3.91113314230973) q[2];
sx q[2];
rz(9.65658632516071) q[2];
rz(1.21094632148743) q[3];
sx q[3];
rz(4.32391932805116) q[3];
sx q[3];
rz(9.7501615345399) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.668321907520294) q[0];
sx q[0];
rz(3.54596916039521) q[0];
sx q[0];
rz(9.80165082811519) q[0];
rz(0.511415243148804) q[1];
sx q[1];
rz(1.8835952599817) q[1];
sx q[1];
rz(11.7223920583646) q[1];
cx q[1],q[0];
rz(0.91746860742569) q[0];
sx q[0];
rz(3.25278029044206) q[0];
sx q[0];
rz(9.58263332246944) q[0];
rz(1.78819441795349) q[2];
sx q[2];
rz(4.1442210992151) q[2];
sx q[2];
rz(11.3142266035001) q[2];
cx q[2],q[1];
rz(-0.626292645931244) q[1];
sx q[1];
rz(3.90230700572068) q[1];
sx q[1];
rz(10.2619667410771) q[1];
rz(0.021447354927659) q[3];
sx q[3];
rz(4.19081583817536) q[3];
sx q[3];
rz(10.0790438413541) q[3];
cx q[3],q[2];
rz(1.97229981422424) q[2];
sx q[2];
rz(4.15631392796571) q[2];
sx q[2];
rz(8.72673729657336) q[2];
rz(1.19971477985382) q[3];
sx q[3];
rz(4.04727658827836) q[3];
sx q[3];
rz(10.0244411587636) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.190716281533241) q[0];
sx q[0];
rz(3.41601279576356) q[0];
sx q[0];
rz(9.55050387083694) q[0];
rz(1.02672863006592) q[1];
sx q[1];
rz(4.89602878888185) q[1];
sx q[1];
rz(9.95098987816974) q[1];
cx q[1],q[0];
rz(0.95740282535553) q[0];
sx q[0];
rz(4.23482778866822) q[0];
sx q[0];
rz(9.34123169480964) q[0];
rz(-0.194296464323997) q[2];
sx q[2];
rz(3.84777233202989) q[2];
sx q[2];
rz(9.88225973247691) q[2];
cx q[2],q[1];
rz(0.806128203868866) q[1];
sx q[1];
rz(4.60996177990968) q[1];
sx q[1];
rz(10.6138753652494) q[1];
rz(0.775330603122711) q[3];
sx q[3];
rz(1.03639999230439) q[3];
sx q[3];
rz(8.86693952082797) q[3];
cx q[3],q[2];
rz(-1.04782629013062) q[2];
sx q[2];
rz(3.97950038512284) q[2];
sx q[2];
rz(9.114518916599) q[2];
rz(0.0349164977669716) q[3];
sx q[3];
rz(4.89713767369325) q[3];
sx q[3];
rz(10.2144807338636) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.07636153697968) q[0];
sx q[0];
rz(4.6236379464441) q[0];
sx q[0];
rz(9.05374798773929) q[0];
rz(0.0648443549871445) q[1];
sx q[1];
rz(3.96386715968186) q[1];
sx q[1];
rz(11.2993684768598) q[1];
cx q[1],q[0];
rz(0.933439135551453) q[0];
sx q[0];
rz(4.36981347401673) q[0];
sx q[0];
rz(10.10571793317) q[0];
rz(1.92239022254944) q[2];
sx q[2];
rz(4.11291977961595) q[2];
sx q[2];
rz(8.77770397662326) q[2];
cx q[2],q[1];
rz(0.632584095001221) q[1];
sx q[1];
rz(1.96880892117555) q[1];
sx q[1];
rz(9.38741311653658) q[1];
rz(1.10478496551514) q[3];
sx q[3];
rz(1.31068411667878) q[3];
sx q[3];
rz(9.08787605761691) q[3];
cx q[3],q[2];
rz(0.123966991901398) q[2];
sx q[2];
rz(4.61831727822358) q[2];
sx q[2];
rz(10.963122701637) q[2];
rz(-0.416878372430801) q[3];
sx q[3];
rz(3.64850518305833) q[3];
sx q[3];
rz(9.59513725935622) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0403195470571518) q[0];
sx q[0];
rz(2.90770010848577) q[0];
sx q[0];
rz(9.85341373681232) q[0];
rz(-1.12426948547363) q[1];
sx q[1];
rz(4.74399569829042) q[1];
sx q[1];
rz(10.1943557619969) q[1];
cx q[1],q[0];
rz(1.29630994796753) q[0];
sx q[0];
rz(3.14519506751607) q[0];
sx q[0];
rz(8.88930062054797) q[0];
rz(-0.495411843061447) q[2];
sx q[2];
rz(4.24253896077211) q[2];
sx q[2];
rz(9.52227145283624) q[2];
cx q[2],q[1];
rz(0.102485485374928) q[1];
sx q[1];
rz(4.41884067853028) q[1];
sx q[1];
rz(10.3744462490003) q[1];
rz(0.510975241661072) q[3];
sx q[3];
rz(1.87124124367768) q[3];
sx q[3];
rz(10.4425580263059) q[3];
cx q[3],q[2];
rz(0.65772420167923) q[2];
sx q[2];
rz(1.85720733006532) q[2];
sx q[2];
rz(9.5202731475155) q[2];
rz(0.541984498500824) q[3];
sx q[3];
rz(2.67545500596101) q[3];
sx q[3];
rz(9.29556148349448) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.995904564857483) q[0];
sx q[0];
rz(4.07197693188722) q[0];
sx q[0];
rz(9.53659666924878) q[0];
rz(1.31891405582428) q[1];
sx q[1];
rz(5.03837886651094) q[1];
sx q[1];
rz(10.730439400665) q[1];
cx q[1],q[0];
rz(2.99667167663574) q[0];
sx q[0];
rz(4.83614364464814) q[0];
sx q[0];
rz(9.97314558028384) q[0];
rz(-1.1131067276001) q[2];
sx q[2];
rz(2.82972422440583) q[2];
sx q[2];
rz(11.517251944534) q[2];
cx q[2],q[1];
rz(-0.240519240498543) q[1];
sx q[1];
rz(1.85226336319978) q[1];
sx q[1];
rz(10.6063656568448) q[1];
rz(0.420283854007721) q[3];
sx q[3];
rz(4.20434919198091) q[3];
sx q[3];
rz(10.2884211897771) q[3];
cx q[3],q[2];
rz(-0.328851133584976) q[2];
sx q[2];
rz(4.29435828526551) q[2];
sx q[2];
rz(10.5444934129636) q[2];
rz(0.268519222736359) q[3];
sx q[3];
rz(4.8790067752176) q[3];
sx q[3];
rz(10.6102895498197) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.280652284622192) q[0];
sx q[0];
rz(3.90345427592332) q[0];
sx q[0];
rz(8.81583408116504) q[0];
rz(0.877434194087982) q[1];
sx q[1];
rz(4.14336326916749) q[1];
sx q[1];
rz(10.0370216131131) q[1];
cx q[1],q[0];
rz(-0.378924816846848) q[0];
sx q[0];
rz(3.98578974803025) q[0];
sx q[0];
rz(8.40569660662814) q[0];
rz(-1.04184091091156) q[2];
sx q[2];
rz(3.77966973383958) q[2];
sx q[2];
rz(11.8660468816678) q[2];
cx q[2],q[1];
rz(0.228371992707253) q[1];
sx q[1];
rz(2.29715839226777) q[1];
sx q[1];
rz(10.3642579674642) q[1];
rz(1.02965712547302) q[3];
sx q[3];
rz(4.64553931553895) q[3];
sx q[3];
rz(9.98614201544925) q[3];
cx q[3],q[2];
rz(-0.0732542425394058) q[2];
sx q[2];
rz(4.74275949795777) q[2];
sx q[2];
rz(9.2373181193988) q[2];
rz(1.11202621459961) q[3];
sx q[3];
rz(4.05968913634355) q[3];
sx q[3];
rz(9.90824360250636) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.407923132181168) q[0];
sx q[0];
rz(3.9556281884485) q[0];
sx q[0];
rz(9.66878305970832) q[0];
rz(1.65812742710114) q[1];
sx q[1];
rz(4.76419857342774) q[1];
sx q[1];
rz(10.1974528193395) q[1];
cx q[1],q[0];
rz(1.34983944892883) q[0];
sx q[0];
rz(2.04492214520509) q[0];
sx q[0];
rz(9.98743090628787) q[0];
rz(1.24638903141022) q[2];
sx q[2];
rz(4.58232370217378) q[2];
sx q[2];
rz(9.03217468260928) q[2];
cx q[2],q[1];
rz(0.150688424706459) q[1];
sx q[1];
rz(3.48782345850999) q[1];
sx q[1];
rz(11.3266383171003) q[1];
rz(0.783873856067657) q[3];
sx q[3];
rz(4.0996944626146) q[3];
sx q[3];
rz(10.1161014199178) q[3];
cx q[3],q[2];
rz(1.33857583999634) q[2];
sx q[2];
rz(4.10279122193391) q[2];
sx q[2];
rz(10.1516639947812) q[2];
rz(1.10729849338531) q[3];
sx q[3];
rz(3.65032699902589) q[3];
sx q[3];
rz(8.41748151778384) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.9299328327179) q[0];
sx q[0];
rz(3.58458924491937) q[0];
sx q[0];
rz(9.03274980782672) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.68397343158722) q[1];
sx q[1];
rz(4.1576431115442) q[1];
sx q[1];
rz(10.0056122898976) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.588883697986603) q[2];
sx q[2];
rz(3.58341852028901) q[2];
sx q[2];
rz(10.5707119464795) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.398409396409988) q[3];
sx q[3];
rz(4.51475134690339) q[3];
sx q[3];
rz(9.56221442519828) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
