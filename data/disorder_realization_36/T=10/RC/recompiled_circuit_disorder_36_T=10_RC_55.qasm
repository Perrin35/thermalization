OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.261895) q[0];
sx q[0];
rz(-1.7356153) q[0];
sx q[0];
rz(0.66189712) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53519997) q[2];
sx q[2];
rz(-2.0173965) q[2];
sx q[2];
rz(-1.610178) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56596245) q[1];
sx q[1];
rz(-1.168025) q[1];
sx q[1];
rz(2.8494542) q[1];
x q[2];
rz(2.203381) q[3];
sx q[3];
rz(-1.0217561) q[3];
sx q[3];
rz(3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082367912) q[0];
sx q[0];
rz(-1.1384283) q[0];
sx q[0];
rz(0.57605497) q[0];
rz(-pi) q[1];
rz(1.4405865) q[2];
sx q[2];
rz(-1.7034334) q[2];
sx q[2];
rz(2.8108033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80026514) q[1];
sx q[1];
rz(-2.6086573) q[1];
sx q[1];
rz(2.76782) q[1];
x q[2];
rz(-0.91652292) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-0.4804002) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-1.0864331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.438293) q[0];
sx q[0];
rz(-2.6801531) q[0];
sx q[0];
rz(-1.567054) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11523192) q[2];
sx q[2];
rz(-2.2149202) q[2];
sx q[2];
rz(0.53172058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4194581) q[1];
sx q[1];
rz(-0.28885435) q[1];
sx q[1];
rz(3.1232749) q[1];
rz(-pi) q[2];
rz(-1.1658737) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(-2.7584372) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.6960467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955129) q[0];
sx q[0];
rz(-2.4582986) q[0];
sx q[0];
rz(-0.28738316) q[0];
rz(1.4807329) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(0.71722523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7062263) q[1];
sx q[1];
rz(-1.5063783) q[1];
sx q[1];
rz(0.21965841) q[1];
rz(-2.0914145) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-3.1307427) q[0];
rz(-pi) q[1];
rz(-0.48798497) q[2];
sx q[2];
rz(-0.93163604) q[2];
sx q[2];
rz(2.0069063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2367868) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(-2.3805815) q[1];
rz(-1.0075931) q[3];
sx q[3];
rz(-2.0493205) q[3];
sx q[3];
rz(-2.0625045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-0.12621005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70996767) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(-2.9931195) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(-2.0691878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0969442) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(-2.2915927) q[1];
x q[2];
rz(1.24228) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643499) q[0];
sx q[0];
rz(-0.61315216) q[0];
sx q[0];
rz(-0.22293798) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26142188) q[2];
sx q[2];
rz(-1.6209507) q[2];
sx q[2];
rz(-2.9068771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4899788) q[1];
sx q[1];
rz(-1.4414806) q[1];
sx q[1];
rz(3.053385) q[1];
rz(-pi) q[2];
rz(1.3748752) q[3];
sx q[3];
rz(-2.19176) q[3];
sx q[3];
rz(1.2697112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-2.7977978) q[2];
rz(2.5750459) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(-2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-0.87160814) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5430785) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(-3.024858) q[0];
rz(-pi) q[1];
rz(-0.34198728) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(-2.4004186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1069113) q[1];
sx q[1];
rz(-0.92202631) q[1];
sx q[1];
rz(-2.0295063) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.4551216) q[3];
sx q[3];
rz(2.065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(0.16601673) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.6171914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6451384) q[0];
sx q[0];
rz(-1.5418515) q[0];
sx q[0];
rz(-1.4883947) q[0];
rz(-pi) q[1];
rz(-2.6644601) q[2];
sx q[2];
rz(-1.4002561) q[2];
sx q[2];
rz(-1.0605304) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9164239) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(1.9041063) q[1];
rz(-pi) q[2];
rz(2.7306261) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(-2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8966658) q[0];
sx q[0];
rz(-1.1765119) q[0];
sx q[0];
rz(3.0299597) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16742736) q[2];
sx q[2];
rz(-1.9055467) q[2];
sx q[2];
rz(-1.6850922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0512143) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(-1.315209) q[1];
x q[2];
rz(-0.6110544) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(2.4886481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(0.28094963) q[2];
sx q[2];
rz(-1.8200257) q[2];
sx q[2];
rz(-0.65489468) q[2];
rz(0.67646277) q[3];
sx q[3];
rz(-0.87920311) q[3];
sx q[3];
rz(1.9999947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
