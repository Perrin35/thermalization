OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(-1.1298236) q[0];
sx q[0];
rz(3.1273754) q[0];
rz(-3.0832503) q[1];
sx q[1];
rz(-0.74906936) q[1];
sx q[1];
rz(-0.60325375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8612807) q[0];
sx q[0];
rz(-2.2076108) q[0];
sx q[0];
rz(-2.0122819) q[0];
x q[1];
rz(1.649734) q[2];
sx q[2];
rz(-0.65627718) q[2];
sx q[2];
rz(-0.62103926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4976884) q[1];
sx q[1];
rz(-0.7086646) q[1];
sx q[1];
rz(-2.9008615) q[1];
rz(-pi) q[2];
rz(-2.2774062) q[3];
sx q[3];
rz(-1.7134616) q[3];
sx q[3];
rz(-0.29970009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76838905) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(1.0962037) q[2];
rz(2.0632035) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(-2.6064742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815149) q[0];
sx q[0];
rz(-2.09477) q[0];
sx q[0];
rz(-0.4775508) q[0];
rz(1.011147) q[1];
sx q[1];
rz(-0.54713455) q[1];
sx q[1];
rz(2.0571041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3318429) q[0];
sx q[0];
rz(-1.4536152) q[0];
sx q[0];
rz(1.940889) q[0];
rz(-pi) q[1];
rz(-1.162503) q[2];
sx q[2];
rz(-1.5973685) q[2];
sx q[2];
rz(-1.0607189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62821416) q[1];
sx q[1];
rz(-2.5625554) q[1];
sx q[1];
rz(-2.8138334) q[1];
rz(-0.22900692) q[3];
sx q[3];
rz(-2.3888616) q[3];
sx q[3];
rz(1.5040042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61887211) q[2];
sx q[2];
rz(-0.42417696) q[2];
sx q[2];
rz(1.0901394) q[2];
rz(-2.5777396) q[3];
sx q[3];
rz(-2.4520051) q[3];
sx q[3];
rz(0.032698154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76145935) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(-1.5367966) q[0];
rz(1.3179294) q[1];
sx q[1];
rz(-2.047796) q[1];
sx q[1];
rz(-2.7440548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1434162) q[0];
sx q[0];
rz(-0.93188028) q[0];
sx q[0];
rz(-2.8692416) q[0];
rz(-pi) q[1];
rz(-0.2961646) q[2];
sx q[2];
rz(-1.0131179) q[2];
sx q[2];
rz(-2.9147343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6315115) q[1];
sx q[1];
rz(-2.1891174) q[1];
sx q[1];
rz(-2.8272998) q[1];
rz(-pi) q[2];
rz(0.52513772) q[3];
sx q[3];
rz(-2.0179837) q[3];
sx q[3];
rz(0.95765985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2354108) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(-0.57514352) q[2];
rz(2.4681674) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(1.5967691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49066082) q[0];
sx q[0];
rz(-2.0199825) q[0];
sx q[0];
rz(0.90890539) q[0];
rz(0.017223651) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-3.1294894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4151013) q[0];
sx q[0];
rz(-0.68031497) q[0];
sx q[0];
rz(2.6136287) q[0];
rz(-0.66476314) q[2];
sx q[2];
rz(-2.3741907) q[2];
sx q[2];
rz(0.27945159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3084532) q[1];
sx q[1];
rz(-0.97452449) q[1];
sx q[1];
rz(-0.27383974) q[1];
rz(3.0148325) q[3];
sx q[3];
rz(-1.6157627) q[3];
sx q[3];
rz(0.25925207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.072731344) q[2];
sx q[2];
rz(-2.8673745) q[2];
sx q[2];
rz(-2.7552628) q[2];
rz(2.7522411) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(1.0079591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35578457) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(0.32387787) q[0];
rz(0.38662275) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(-1.5868384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64774871) q[0];
sx q[0];
rz(-0.84058981) q[0];
sx q[0];
rz(1.4842008) q[0];
rz(-pi) q[1];
rz(1.0350758) q[2];
sx q[2];
rz(-1.8958099) q[2];
sx q[2];
rz(2.6485788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91134478) q[1];
sx q[1];
rz(-0.94083386) q[1];
sx q[1];
rz(-2.1542284) q[1];
rz(0.50916785) q[3];
sx q[3];
rz(-0.8664136) q[3];
sx q[3];
rz(2.1455163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7524449) q[2];
sx q[2];
rz(-0.84263313) q[2];
sx q[2];
rz(-3.0774934) q[2];
rz(-2.5373503) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(1.2324246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-2.5305643) q[0];
sx q[0];
rz(0.87919277) q[0];
rz(2.7912256) q[1];
sx q[1];
rz(-1.3060952) q[1];
sx q[1];
rz(-0.15215692) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76497117) q[0];
sx q[0];
rz(-0.76036727) q[0];
sx q[0];
rz(2.7432824) q[0];
rz(-1.9273571) q[2];
sx q[2];
rz(-2.6313553) q[2];
sx q[2];
rz(1.4067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81318839) q[1];
sx q[1];
rz(-2.0021353) q[1];
sx q[1];
rz(-1.793641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.426126) q[3];
sx q[3];
rz(-2.5757534) q[3];
sx q[3];
rz(-2.6239606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0396314) q[2];
sx q[2];
rz(-2.8122718) q[2];
sx q[2];
rz(0.24001089) q[2];
rz(0.65252423) q[3];
sx q[3];
rz(-2.0034761) q[3];
sx q[3];
rz(1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0535102) q[0];
sx q[0];
rz(-1.9048012) q[0];
sx q[0];
rz(0.85813338) q[0];
rz(1.3872604) q[1];
sx q[1];
rz(-2.6871082) q[1];
sx q[1];
rz(0.33285704) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4384321) q[0];
sx q[0];
rz(-2.1370208) q[0];
sx q[0];
rz(-2.7661408) q[0];
rz(1.5120686) q[2];
sx q[2];
rz(-1.3212122) q[2];
sx q[2];
rz(-1.5098315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.27204) q[1];
sx q[1];
rz(-1.365188) q[1];
sx q[1];
rz(-0.42821347) q[1];
rz(-0.65471411) q[3];
sx q[3];
rz(-1.8548449) q[3];
sx q[3];
rz(-0.14446196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7994999) q[2];
sx q[2];
rz(-0.19793333) q[2];
sx q[2];
rz(-1.0510772) q[2];
rz(1.6445232) q[3];
sx q[3];
rz(-1.9688508) q[3];
sx q[3];
rz(0.79819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368211) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(2.2494702) q[0];
rz(0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(-2.3687252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.597201) q[0];
sx q[0];
rz(-2.3344451) q[0];
sx q[0];
rz(-2.3788484) q[0];
x q[1];
rz(-0.70336999) q[2];
sx q[2];
rz(-1.7711319) q[2];
sx q[2];
rz(-3.1205265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9443431) q[1];
sx q[1];
rz(-2.5641003) q[1];
sx q[1];
rz(0.61402278) q[1];
rz(-0.77247932) q[3];
sx q[3];
rz(-2.8101343) q[3];
sx q[3];
rz(2.9951688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.2246776) q[2];
sx q[2];
rz(2.4490855) q[2];
rz(-2.6912189) q[3];
sx q[3];
rz(-1.2053442) q[3];
sx q[3];
rz(-1.6374121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57167989) q[0];
sx q[0];
rz(-0.79140651) q[0];
sx q[0];
rz(2.1929542) q[0];
rz(2.0750849) q[1];
sx q[1];
rz(-2.152161) q[1];
sx q[1];
rz(-1.271064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1656837) q[0];
sx q[0];
rz(-1.6635487) q[0];
sx q[0];
rz(-1.3646056) q[0];
x q[1];
rz(-0.99540751) q[2];
sx q[2];
rz(-0.040608309) q[2];
sx q[2];
rz(0.81088582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3349325) q[1];
sx q[1];
rz(-2.1927823) q[1];
sx q[1];
rz(0.16652624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4404092) q[3];
sx q[3];
rz(-2.009863) q[3];
sx q[3];
rz(1.8177123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5991685) q[2];
sx q[2];
rz(-1.3072661) q[2];
sx q[2];
rz(-3.001281) q[2];
rz(-3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(-1.2062581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62823826) q[0];
sx q[0];
rz(-1.6980549) q[0];
sx q[0];
rz(-0.7775318) q[0];
rz(0.93742049) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(0.44313988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454058) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(3.0565673) q[0];
x q[1];
rz(-0.32989008) q[2];
sx q[2];
rz(-2.1426845) q[2];
sx q[2];
rz(3.013333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7378583) q[1];
sx q[1];
rz(-2.1843806) q[1];
sx q[1];
rz(-1.1627083) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0983766) q[3];
sx q[3];
rz(-2.7948423) q[3];
sx q[3];
rz(-2.7048793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-0.12760663) q[2];
rz(2.8429032) q[3];
sx q[3];
rz(-1.1726215) q[3];
sx q[3];
rz(-2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5304607) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(2.5313189) q[1];
sx q[1];
rz(-0.88519575) q[1];
sx q[1];
rz(-2.0559678) q[1];
rz(-1.4971785) q[2];
sx q[2];
rz(-0.72009077) q[2];
sx q[2];
rz(2.1511267) q[2];
rz(-2.0213303) q[3];
sx q[3];
rz(-1.3792752) q[3];
sx q[3];
rz(-3.006912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
