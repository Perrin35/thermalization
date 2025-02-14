OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38464889) q[0];
sx q[0];
rz(-2.3031213) q[0];
sx q[0];
rz(-2.5323001) q[0];
rz(-2.4203909) q[1];
sx q[1];
rz(-2.2130122) q[1];
sx q[1];
rz(2.0452926) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247835) q[0];
sx q[0];
rz(-1.4583298) q[0];
sx q[0];
rz(-0.87392585) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47639552) q[2];
sx q[2];
rz(-0.98586776) q[2];
sx q[2];
rz(-0.75854075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9043353) q[1];
sx q[1];
rz(-1.0544027) q[1];
sx q[1];
rz(1.2468491) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5240229) q[3];
sx q[3];
rz(-2.9520116) q[3];
sx q[3];
rz(-1.1763587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91813749) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(-1.9523331) q[2];
rz(0.39092815) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(1.1322359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.3005715) q[1];
sx q[1];
rz(0.71151412) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92828099) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(0.9071945) q[0];
rz(-pi) q[1];
rz(2.6911084) q[2];
sx q[2];
rz(-0.50651359) q[2];
sx q[2];
rz(-1.5114558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1700279) q[1];
sx q[1];
rz(-1.6508647) q[1];
sx q[1];
rz(-0.34645924) q[1];
x q[2];
rz(0.94717104) q[3];
sx q[3];
rz(-1.5963364) q[3];
sx q[3];
rz(0.21835777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(-2.0785275) q[2];
rz(1.6478018) q[3];
sx q[3];
rz(-0.46896514) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4110334) q[0];
sx q[0];
rz(-2.7506802) q[0];
sx q[0];
rz(-2.539769) q[0];
rz(2.9959294) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(0.62785968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74341853) q[0];
sx q[0];
rz(-1.6891915) q[0];
sx q[0];
rz(1.5009686) q[0];
rz(-2.7854087) q[2];
sx q[2];
rz(-2.0938452) q[2];
sx q[2];
rz(1.5580391) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9380629) q[1];
sx q[1];
rz(-0.4258241) q[1];
sx q[1];
rz(0.99217207) q[1];
rz(-pi) q[2];
rz(-0.95471621) q[3];
sx q[3];
rz(-2.4297415) q[3];
sx q[3];
rz(2.4095132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38498983) q[2];
sx q[2];
rz(-1.0641229) q[2];
sx q[2];
rz(-0.21558726) q[2];
rz(-1.138341) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1348006) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(-3.0227645) q[0];
rz(1.6670594) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(0.80783358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87619239) q[0];
sx q[0];
rz(-1.0256488) q[0];
sx q[0];
rz(-0.67727526) q[0];
rz(-pi) q[1];
rz(-3.0962871) q[2];
sx q[2];
rz(-1.2477562) q[2];
sx q[2];
rz(-1.3458061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68271676) q[1];
sx q[1];
rz(-1.7380875) q[1];
sx q[1];
rz(0.95750995) q[1];
rz(3.1224072) q[3];
sx q[3];
rz(-0.67836715) q[3];
sx q[3];
rz(0.32103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5011751) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(-0.51587063) q[2];
rz(3.0473895) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(-1.5883821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.9544253) q[0];
rz(0.59133235) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(2.3147413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39385228) q[0];
sx q[0];
rz(-1.4497787) q[0];
sx q[0];
rz(0.38328247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0127147) q[2];
sx q[2];
rz(-1.5019755) q[2];
sx q[2];
rz(-2.0134913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34870711) q[1];
sx q[1];
rz(-1.4408083) q[1];
sx q[1];
rz(-2.2568364) q[1];
rz(-pi) q[2];
rz(0.99991344) q[3];
sx q[3];
rz(-1.1383082) q[3];
sx q[3];
rz(1.692358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4225215) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(-2.1264326) q[2];
rz(0.90967956) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5231617) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(0.30995187) q[0];
rz(3.1228206) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(0.087336691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76742889) q[0];
sx q[0];
rz(-1.5383676) q[0];
sx q[0];
rz(-0.57291605) q[0];
rz(2.3051065) q[2];
sx q[2];
rz(-2.2936471) q[2];
sx q[2];
rz(-1.370887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0897905) q[1];
sx q[1];
rz(-0.6615122) q[1];
sx q[1];
rz(-2.2480985) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9841626) q[3];
sx q[3];
rz(-1.0478019) q[3];
sx q[3];
rz(-0.98539814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(3.0067048) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(-2.7556162) q[0];
rz(2.0416253) q[1];
sx q[1];
rz(-0.67953449) q[1];
sx q[1];
rz(0.78752548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.814491) q[0];
sx q[0];
rz(-1.2145894) q[0];
sx q[0];
rz(2.8339631) q[0];
rz(-pi) q[1];
rz(1.196127) q[2];
sx q[2];
rz(-2.154899) q[2];
sx q[2];
rz(2.0541592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0037076) q[1];
sx q[1];
rz(-0.862993) q[1];
sx q[1];
rz(-2.6801609) q[1];
rz(1.8058895) q[3];
sx q[3];
rz(-1.4054148) q[3];
sx q[3];
rz(2.7602285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(1.5137399) q[2];
rz(-1.2210023) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(-0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8319703) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(0.23304644) q[0];
rz(-1.6538992) q[1];
sx q[1];
rz(-2.1079886) q[1];
sx q[1];
rz(-2.1344562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3171948) q[0];
sx q[0];
rz(-1.2318512) q[0];
sx q[0];
rz(2.5934153) q[0];
x q[1];
rz(0.39733823) q[2];
sx q[2];
rz(-0.66911829) q[2];
sx q[2];
rz(-0.64815176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1888652) q[1];
sx q[1];
rz(-0.42546526) q[1];
sx q[1];
rz(-1.9317606) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5910025) q[3];
sx q[3];
rz(-1.8148225) q[3];
sx q[3];
rz(3.1283825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2795589) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(-0.88279185) q[2];
rz(-0.27859303) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(-2.6826503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19926628) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(-1.5698154) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-2.4615361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5424283) q[0];
sx q[0];
rz(-2.3788646) q[0];
sx q[0];
rz(3.0551856) q[0];
rz(-1.6974738) q[2];
sx q[2];
rz(-1.6886204) q[2];
sx q[2];
rz(0.5165259) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16526651) q[1];
sx q[1];
rz(-1.191947) q[1];
sx q[1];
rz(0.59013175) q[1];
rz(0.84362883) q[3];
sx q[3];
rz(-1.9988201) q[3];
sx q[3];
rz(2.6035863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.9752768) q[2];
rz(1.204528) q[3];
sx q[3];
rz(-1.8185936) q[3];
sx q[3];
rz(-2.5148463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.523943) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(-0.43706617) q[0];
rz(2.3433459) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(-2.9214568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6740538) q[0];
sx q[0];
rz(-2.3706145) q[0];
sx q[0];
rz(1.6378372) q[0];
rz(-pi) q[1];
x q[1];
rz(1.885845) q[2];
sx q[2];
rz(-1.7654668) q[2];
sx q[2];
rz(0.48182975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70500532) q[1];
sx q[1];
rz(-1.2975177) q[1];
sx q[1];
rz(-2.336034) q[1];
rz(3.0978904) q[3];
sx q[3];
rz(-0.80207295) q[3];
sx q[3];
rz(-1.613783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2685214) q[2];
sx q[2];
rz(-1.9067418) q[2];
sx q[2];
rz(-2.851167) q[2];
rz(0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.2344454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3270522) q[0];
sx q[0];
rz(-0.70801133) q[0];
sx q[0];
rz(1.1109362) q[0];
rz(0.064432714) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(0.031485166) q[2];
sx q[2];
rz(-0.97618547) q[2];
sx q[2];
rz(0.57913274) q[2];
rz(1.8364344) q[3];
sx q[3];
rz(-2.0943741) q[3];
sx q[3];
rz(1.5423273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
