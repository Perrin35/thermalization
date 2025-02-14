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
rz(-3.0120612) q[0];
sx q[0];
rz(-3.01053) q[0];
sx q[0];
rz(0.60453209) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(2.8837535) q[1];
sx q[1];
rz(16.937994) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2211232) q[0];
sx q[0];
rz(-1.2386432) q[0];
sx q[0];
rz(1.3060547) q[0];
x q[1];
rz(-2.8240439) q[2];
sx q[2];
rz(-0.63830909) q[2];
sx q[2];
rz(0.93250634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9166419) q[1];
sx q[1];
rz(-1.0826546) q[1];
sx q[1];
rz(-2.3596838) q[1];
x q[2];
rz(-1.3377519) q[3];
sx q[3];
rz(-2.6342158) q[3];
sx q[3];
rz(-0.22358596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0543542) q[2];
sx q[2];
rz(-1.9638991) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-0.68322244) q[3];
sx q[3];
rz(0.78342485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0403274) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-0.26327565) q[0];
rz(-2.1030262) q[1];
sx q[1];
rz(-0.34663215) q[1];
sx q[1];
rz(-2.3668049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84768554) q[0];
sx q[0];
rz(-2.4890702) q[0];
sx q[0];
rz(-2.7987665) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69157522) q[2];
sx q[2];
rz(-2.1082768) q[2];
sx q[2];
rz(-0.36851685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6331433) q[1];
sx q[1];
rz(-1.5337428) q[1];
sx q[1];
rz(-1.9502742) q[1];
x q[2];
rz(-0.56258308) q[3];
sx q[3];
rz(-0.78662164) q[3];
sx q[3];
rz(-2.7275769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7210377) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(2.5431385) q[2];
rz(-1.7789486) q[3];
sx q[3];
rz(-2.2612031) q[3];
sx q[3];
rz(-0.27073282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28629985) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(1.5077952) q[0];
rz(-0.15448054) q[1];
sx q[1];
rz(-1.4198317) q[1];
sx q[1];
rz(2.3522164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76217795) q[0];
sx q[0];
rz(-1.6573849) q[0];
sx q[0];
rz(-0.53126727) q[0];
rz(-0.28114281) q[2];
sx q[2];
rz(-1.4340377) q[2];
sx q[2];
rz(0.006055486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4142609) q[1];
sx q[1];
rz(-1.2664898) q[1];
sx q[1];
rz(-0.93549606) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79929914) q[3];
sx q[3];
rz(-1.2033312) q[3];
sx q[3];
rz(-0.41900837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(1.6645128) q[2];
rz(1.2137671) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(-1.414813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-2.0857683) q[0];
sx q[0];
rz(0.61071998) q[0];
rz(-2.4702813) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(-1.3549365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5089534) q[0];
sx q[0];
rz(-1.1221948) q[0];
sx q[0];
rz(-0.53855702) q[0];
rz(-pi) q[1];
rz(-1.1995537) q[2];
sx q[2];
rz(-1.4535731) q[2];
sx q[2];
rz(-3.0407112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0555041) q[1];
sx q[1];
rz(-1.633606) q[1];
sx q[1];
rz(1.8163866) q[1];
x q[2];
rz(1.7132961) q[3];
sx q[3];
rz(-0.5915407) q[3];
sx q[3];
rz(1.174389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2063724) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(-0.26885232) q[2];
rz(-0.74784652) q[3];
sx q[3];
rz(-2.3220389) q[3];
sx q[3];
rz(1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.185323) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(-2.4638033) q[0];
rz(-0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.6544624) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8003755) q[0];
sx q[0];
rz(-1.7066597) q[0];
sx q[0];
rz(-0.02639833) q[0];
rz(-pi) q[1];
rz(1.1688091) q[2];
sx q[2];
rz(-1.6239407) q[2];
sx q[2];
rz(-2.1316018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9095535) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(-2.1559245) q[1];
x q[2];
rz(-0.048154496) q[3];
sx q[3];
rz(-1.1954525) q[3];
sx q[3];
rz(-1.0182235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80374485) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(-0.35550508) q[2];
rz(1.0454987) q[3];
sx q[3];
rz(-1.9020566) q[3];
sx q[3];
rz(-0.56226468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1657555) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(-0.088706644) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.8314223) q[1];
sx q[1];
rz(0.075693695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3189118) q[0];
sx q[0];
rz(-1.7633798) q[0];
sx q[0];
rz(-0.67355021) q[0];
rz(-1.6275109) q[2];
sx q[2];
rz(-2.4091588) q[2];
sx q[2];
rz(-0.54679856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10403216) q[1];
sx q[1];
rz(-2.1227897) q[1];
sx q[1];
rz(-2.9747573) q[1];
x q[2];
rz(-1.2457341) q[3];
sx q[3];
rz(-2.1233798) q[3];
sx q[3];
rz(2.0393275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82464108) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(-2.1023882) q[2];
rz(2.9292987) q[3];
sx q[3];
rz(-1.1024691) q[3];
sx q[3];
rz(1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067326389) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(-1.0443895) q[0];
rz(-2.3692865) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(2.4078802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4348616) q[0];
sx q[0];
rz(-1.0821663) q[0];
sx q[0];
rz(0.74669331) q[0];
rz(1.9416503) q[2];
sx q[2];
rz(-0.39526734) q[2];
sx q[2];
rz(-0.42759174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.305334) q[1];
sx q[1];
rz(-1.6871258) q[1];
sx q[1];
rz(-0.29246026) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83603199) q[3];
sx q[3];
rz(-1.6510909) q[3];
sx q[3];
rz(-1.7952331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5250728) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(2.5329242) q[2];
rz(2.0742553) q[3];
sx q[3];
rz(-0.87456861) q[3];
sx q[3];
rz(-0.55060351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50255018) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.4962366) q[0];
rz(1.3642338) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-2.7552557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7832121) q[0];
sx q[0];
rz(-2.2370501) q[0];
sx q[0];
rz(1.9080286) q[0];
x q[1];
rz(1.213721) q[2];
sx q[2];
rz(-0.197796) q[2];
sx q[2];
rz(-1.3551499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5959634) q[1];
sx q[1];
rz(-1.9013604) q[1];
sx q[1];
rz(1.7002374) q[1];
x q[2];
rz(1.6101077) q[3];
sx q[3];
rz(-0.98061258) q[3];
sx q[3];
rz(-3.085768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51314696) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(2.9885542) q[2];
rz(-0.89093527) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(-0.74131596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1560169) q[0];
sx q[0];
rz(-2.9122536) q[0];
sx q[0];
rz(-2.7942221) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(-0.52245021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.98451) q[0];
sx q[0];
rz(-1.9831796) q[0];
sx q[0];
rz(0.097481485) q[0];
rz(1.5711196) q[2];
sx q[2];
rz(-1.9362893) q[2];
sx q[2];
rz(1.5207639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3047239) q[1];
sx q[1];
rz(-0.68075276) q[1];
sx q[1];
rz(-0.70226837) q[1];
rz(-pi) q[2];
rz(-2.7344875) q[3];
sx q[3];
rz(-2.0927755) q[3];
sx q[3];
rz(2.2996421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5106421) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(-1.2003215) q[2];
rz(-0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(-0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2931622) q[0];
sx q[0];
rz(-2.1577305) q[0];
sx q[0];
rz(1.7278607) q[0];
rz(3.0005455) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(-0.65473762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4439693) q[0];
sx q[0];
rz(-1.807741) q[0];
sx q[0];
rz(2.9276642) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0789924) q[2];
sx q[2];
rz(-2.104665) q[2];
sx q[2];
rz(1.2956308) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84147553) q[1];
sx q[1];
rz(-1.3285331) q[1];
sx q[1];
rz(0.17657847) q[1];
rz(-pi) q[2];
x q[2];
rz(1.039417) q[3];
sx q[3];
rz(-1.0194885) q[3];
sx q[3];
rz(-1.2599864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.1620136) q[2];
sx q[2];
rz(-2.4164825) q[2];
rz(-0.890598) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643322) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(2.019885) q[1];
sx q[1];
rz(-1.5678761) q[1];
sx q[1];
rz(-0.97490464) q[1];
rz(2.4779392) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(1.411047) q[3];
sx q[3];
rz(-0.86426576) q[3];
sx q[3];
rz(-0.15440253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
