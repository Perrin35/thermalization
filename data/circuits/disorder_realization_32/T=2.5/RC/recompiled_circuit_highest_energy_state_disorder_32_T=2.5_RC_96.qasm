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
rz(0.12953144) q[0];
sx q[0];
rz(-0.13106267) q[0];
sx q[0];
rz(-0.60453209) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(-0.25783917) q[1];
sx q[1];
rz(-1.2300307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9204694) q[0];
sx q[0];
rz(-1.9029494) q[0];
sx q[0];
rz(-1.3060547) q[0];
rz(-pi) q[1];
rz(-2.8240439) q[2];
sx q[2];
rz(-2.5032836) q[2];
sx q[2];
rz(2.2090863) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22495077) q[1];
sx q[1];
rz(-2.058938) q[1];
sx q[1];
rz(-2.3596838) q[1];
rz(-1.0749726) q[3];
sx q[3];
rz(-1.6832441) q[3];
sx q[3];
rz(1.551764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0872385) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(-3.0053906) q[2];
rz(-0.44998351) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1012652) q[0];
sx q[0];
rz(-0.7985006) q[0];
sx q[0];
rz(2.878317) q[0];
rz(1.0385665) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(2.3668049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2939071) q[0];
sx q[0];
rz(-2.4890702) q[0];
sx q[0];
rz(2.7987665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4500174) q[2];
sx q[2];
rz(-1.0333158) q[2];
sx q[2];
rz(0.36851685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6331433) q[1];
sx q[1];
rz(-1.5337428) q[1];
sx q[1];
rz(-1.9502742) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5790096) q[3];
sx q[3];
rz(-0.78662164) q[3];
sx q[3];
rz(0.41401573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.420555) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(2.5431385) q[2];
rz(1.362644) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(-2.8708598) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(1.6337974) q[0];
rz(2.9871121) q[1];
sx q[1];
rz(-1.4198317) q[1];
sx q[1];
rz(2.3522164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75784439) q[0];
sx q[0];
rz(-1.0417291) q[0];
sx q[0];
rz(1.6711414) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8604498) q[2];
sx q[2];
rz(-1.707555) q[2];
sx q[2];
rz(3.1355372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7273318) q[1];
sx q[1];
rz(-1.8751029) q[1];
sx q[1];
rz(0.93549606) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0752729) q[3];
sx q[3];
rz(-0.83791997) q[3];
sx q[3];
rz(-1.5056226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(-1.6645128) q[2];
rz(1.9278256) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(1.414813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690935) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(-0.61071998) q[0];
rz(-2.4702813) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(1.7866561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31025654) q[0];
sx q[0];
rz(-2.455157) q[0];
sx q[0];
rz(2.3879334) q[0];
x q[1];
rz(0.12570555) q[2];
sx q[2];
rz(-1.9393688) q[2];
sx q[2];
rz(1.6261795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.902261) q[1];
sx q[1];
rz(-2.888255) q[1];
sx q[1];
rz(-1.8239267) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0464744) q[3];
sx q[3];
rz(-0.98604938) q[3];
sx q[3];
rz(-1.0032391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9352202) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(0.74784652) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95626962) q[0];
sx q[0];
rz(-1.7522426) q[0];
sx q[0];
rz(-0.67778936) q[0];
rz(-0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(-1.4871303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8003755) q[0];
sx q[0];
rz(-1.7066597) q[0];
sx q[0];
rz(0.02639833) q[0];
rz(-pi) q[1];
rz(1.1688091) q[2];
sx q[2];
rz(-1.5176519) q[2];
sx q[2];
rz(2.1316018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2525874) q[1];
sx q[1];
rz(-0.72786268) q[1];
sx q[1];
rz(0.83830203) q[1];
x q[2];
rz(-3.0934382) q[3];
sx q[3];
rz(-1.1954525) q[3];
sx q[3];
rz(-2.1233692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3378478) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(-2.7860876) q[2];
rz(-1.0454987) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(-0.56226468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97583714) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(-0.088706644) q[0];
rz(1.5345796) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(3.065899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189118) q[0];
sx q[0];
rz(-1.3782129) q[0];
sx q[0];
rz(2.4680424) q[0];
rz(1.5140818) q[2];
sx q[2];
rz(-0.73243388) q[2];
sx q[2];
rz(-2.5947941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20670375) q[1];
sx q[1];
rz(-2.5674501) q[1];
sx q[1];
rz(1.3074134) q[1];
rz(2.5646943) q[3];
sx q[3];
rz(-1.8460974) q[3];
sx q[3];
rz(-2.4979765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3169516) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(2.1023882) q[2];
rz(-0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(1.645741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067326389) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(-2.0972032) q[0];
rz(0.77230612) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-0.73371249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334856) q[0];
sx q[0];
rz(-2.2758188) q[0];
sx q[0];
rz(0.66410983) q[0];
rz(-pi) q[1];
rz(1.9416503) q[2];
sx q[2];
rz(-0.39526734) q[2];
sx q[2];
rz(-0.42759174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.305334) q[1];
sx q[1];
rz(-1.4544669) q[1];
sx q[1];
rz(-2.8491324) q[1];
rz(-pi) q[2];
rz(-0.10802631) q[3];
sx q[3];
rz(-2.3026534) q[3];
sx q[3];
rz(2.9895003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61651984) q[2];
sx q[2];
rz(-0.51023444) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(0.38633698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3583806) q[0];
sx q[0];
rz(-0.9045426) q[0];
sx q[0];
rz(-1.9080286) q[0];
rz(1.213721) q[2];
sx q[2];
rz(-2.9437967) q[2];
sx q[2];
rz(-1.7864428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2145591) q[1];
sx q[1];
rz(-2.7874569) q[1];
sx q[1];
rz(-0.35978364) q[1];
x q[2];
rz(0.058606996) q[3];
sx q[3];
rz(-2.5502565) q[3];
sx q[3];
rz(3.1268596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51314696) q[2];
sx q[2];
rz(-1.7690423) q[2];
sx q[2];
rz(2.9885542) q[2];
rz(2.2506574) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(-2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9855758) q[0];
sx q[0];
rz(-2.9122536) q[0];
sx q[0];
rz(2.7942221) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.9719351) q[1];
sx q[1];
rz(-2.6191424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1570826) q[0];
sx q[0];
rz(-1.9831796) q[0];
sx q[0];
rz(-3.0441112) q[0];
rz(-pi) q[1];
rz(-2.7760997) q[2];
sx q[2];
rz(-1.5710982) q[2];
sx q[2];
rz(0.049916849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3047239) q[1];
sx q[1];
rz(-2.4608399) q[1];
sx q[1];
rz(-2.4393243) q[1];
rz(-2.1736828) q[3];
sx q[3];
rz(-0.65015745) q[3];
sx q[3];
rz(-3.0126743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5106421) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(-2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2931622) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(1.4137319) q[0];
rz(-3.0005455) q[1];
sx q[1];
rz(-1.7660564) q[1];
sx q[1];
rz(2.486855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4439693) q[0];
sx q[0];
rz(-1.807741) q[0];
sx q[0];
rz(-2.9276642) q[0];
x q[1];
rz(2.0789924) q[2];
sx q[2];
rz(-2.104665) q[2];
sx q[2];
rz(-1.8459619) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6605986) q[1];
sx q[1];
rz(-0.29876041) q[1];
sx q[1];
rz(0.95282747) q[1];
x q[2];
rz(-2.1021757) q[3];
sx q[3];
rz(-2.1221042) q[3];
sx q[3];
rz(1.2599864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(-2.4164825) q[2];
rz(-2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17726041) q[0];
sx q[0];
rz(-1.6163419) q[0];
sx q[0];
rz(1.1905715) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(-2.6910838) q[2];
sx q[2];
rz(-2.4263739) q[2];
sx q[2];
rz(0.77745773) q[2];
rz(0.71286451) q[3];
sx q[3];
rz(-1.449507) q[3];
sx q[3];
rz(1.5206152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
