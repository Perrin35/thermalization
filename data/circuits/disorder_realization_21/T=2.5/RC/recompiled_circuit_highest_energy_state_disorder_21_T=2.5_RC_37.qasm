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
rz(-1.4323956) q[0];
sx q[0];
rz(-0.30269912) q[0];
sx q[0];
rz(-2.1085289) q[0];
rz(-1.2248224) q[1];
sx q[1];
rz(-1.4900102) q[1];
sx q[1];
rz(2.9472247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9308192) q[0];
sx q[0];
rz(-1.418485) q[0];
sx q[0];
rz(-0.62893008) q[0];
rz(2.7871903) q[2];
sx q[2];
rz(-2.4131219) q[2];
sx q[2];
rz(0.013106339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3857419) q[1];
sx q[1];
rz(-0.057690851) q[1];
sx q[1];
rz(2.2123446) q[1];
x q[2];
rz(1.1367873) q[3];
sx q[3];
rz(-0.99054407) q[3];
sx q[3];
rz(-2.2541719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4091461) q[2];
sx q[2];
rz(-1.1017841) q[2];
sx q[2];
rz(2.5771602) q[2];
rz(-1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.0821575) q[0];
sx q[0];
rz(-0.0075639021) q[0];
sx q[0];
rz(2.2285158) q[0];
rz(2.4572241) q[1];
sx q[1];
rz(-3.1412536) q[1];
sx q[1];
rz(-0.93904644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9371999) q[0];
sx q[0];
rz(-2.4449722) q[0];
sx q[0];
rz(-2.533769) q[0];
rz(-pi) q[1];
rz(-1.5576511) q[2];
sx q[2];
rz(-1.1068543) q[2];
sx q[2];
rz(2.9090403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1582271) q[1];
sx q[1];
rz(-1.905113) q[1];
sx q[1];
rz(2.0086221) q[1];
rz(-1.4165808) q[3];
sx q[3];
rz(-0.41114488) q[3];
sx q[3];
rz(-1.7257041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(-0.89246559) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40176323) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(-2.7318562) q[0];
rz(-3.1306664) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(-1.797537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7342012) q[0];
sx q[0];
rz(-1.0599066) q[0];
sx q[0];
rz(0.023500806) q[0];
rz(-3.0777099) q[2];
sx q[2];
rz(-1.6860355) q[2];
sx q[2];
rz(-1.921953) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5508268) q[1];
sx q[1];
rz(-1.8286341) q[1];
sx q[1];
rz(-0.067036585) q[1];
rz(-pi) q[2];
rz(-1.4555172) q[3];
sx q[3];
rz(-1.56101) q[3];
sx q[3];
rz(-2.9392795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.672013) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(3.0984042) q[2];
rz(1.1399266) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7286872) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(1.8034978) q[0];
rz(2.4463704) q[1];
sx q[1];
rz(-0.1278563) q[1];
sx q[1];
rz(2.7241838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5240751) q[0];
sx q[0];
rz(-1.39912) q[0];
sx q[0];
rz(2.2508932) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2222671) q[2];
sx q[2];
rz(-0.67643316) q[2];
sx q[2];
rz(1.2471022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96046916) q[1];
sx q[1];
rz(-0.6746909) q[1];
sx q[1];
rz(-1.8974278) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2715075) q[3];
sx q[3];
rz(-2.6911754) q[3];
sx q[3];
rz(1.3189486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7030316) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(-1.7523127) q[2];
rz(-2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(-2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97838068) q[0];
sx q[0];
rz(-0.037540171) q[0];
sx q[0];
rz(0.96330825) q[0];
rz(2.9523201) q[1];
sx q[1];
rz(-0.015733868) q[1];
sx q[1];
rz(-2.9761369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2509569) q[0];
sx q[0];
rz(-1.5444618) q[0];
sx q[0];
rz(3.0980397) q[0];
rz(-pi) q[1];
rz(-0.87617434) q[2];
sx q[2];
rz(-1.8037705) q[2];
sx q[2];
rz(2.7829952) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7041118) q[1];
sx q[1];
rz(-2.3299814) q[1];
sx q[1];
rz(-0.1603756) q[1];
x q[2];
rz(-2.8235196) q[3];
sx q[3];
rz(-1.9847365) q[3];
sx q[3];
rz(1.3596168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.307622) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(-1.0597672) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5237913) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(-1.6049438) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(1.4401999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337517) q[0];
sx q[0];
rz(-1.4086723) q[0];
sx q[0];
rz(-2.9723245) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9199706) q[2];
sx q[2];
rz(-1.3563507) q[2];
sx q[2];
rz(-2.1943486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.837484) q[1];
sx q[1];
rz(-0.15481259) q[1];
sx q[1];
rz(0.047708851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0470508) q[3];
sx q[3];
rz(-0.19652995) q[3];
sx q[3];
rz(-2.9158338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3692533) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(0.079027979) q[2];
rz(-2.2760271) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(-2.1178093) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2257776) q[0];
sx q[0];
rz(-3.1362035) q[0];
sx q[0];
rz(-2.916577) q[0];
rz(-2.8354538) q[1];
sx q[1];
rz(-0.016751079) q[1];
sx q[1];
rz(-0.87047815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35246145) q[0];
sx q[0];
rz(-1.5232183) q[0];
sx q[0];
rz(-0.0713047) q[0];
rz(-pi) q[1];
rz(-2.5092431) q[2];
sx q[2];
rz(-0.82201695) q[2];
sx q[2];
rz(2.1123304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5815818) q[1];
sx q[1];
rz(-2.3696179) q[1];
sx q[1];
rz(3.0549269) q[1];
rz(-0.87865307) q[3];
sx q[3];
rz(-0.60135287) q[3];
sx q[3];
rz(0.9200615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6656782) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(2.1062984) q[2];
rz(2.3546442) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(-0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(3.1159478) q[0];
rz(1.8772839) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(2.4425676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8549439) q[0];
sx q[0];
rz(-2.790745) q[0];
sx q[0];
rz(-2.4445266) q[0];
rz(-0.46866663) q[2];
sx q[2];
rz(-1.7231517) q[2];
sx q[2];
rz(-3.0909757) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0009202) q[1];
sx q[1];
rz(-2.7503221) q[1];
sx q[1];
rz(2.9468582) q[1];
rz(2.5949536) q[3];
sx q[3];
rz(-0.99144236) q[3];
sx q[3];
rz(2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8890624) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(2.8177281) q[2];
rz(-2.9845386) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(-2.9878555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(-0.59004849) q[0];
rz(0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-0.86729008) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27408263) q[0];
sx q[0];
rz(-1.6757586) q[0];
sx q[0];
rz(0.63469736) q[0];
rz(-pi) q[1];
rz(-2.1731253) q[2];
sx q[2];
rz(-1.3545389) q[2];
sx q[2];
rz(2.575084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54360753) q[1];
sx q[1];
rz(-1.6389009) q[1];
sx q[1];
rz(1.4582514) q[1];
rz(0.85764472) q[3];
sx q[3];
rz(-1.4271171) q[3];
sx q[3];
rz(-1.194467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1934293) q[2];
sx q[2];
rz(-0.89274222) q[2];
sx q[2];
rz(3.0286922) q[2];
rz(-2.0924977) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(0.73672867) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.076544) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(1.6381868) q[0];
rz(-0.63649559) q[1];
sx q[1];
rz(-2.681585) q[1];
sx q[1];
rz(1.5696625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2777247) q[0];
sx q[0];
rz(-1.6536599) q[0];
sx q[0];
rz(3.0623097) q[0];
x q[1];
rz(-3.1384597) q[2];
sx q[2];
rz(-1.5701446) q[2];
sx q[2];
rz(-2.3359483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0048945) q[1];
sx q[1];
rz(-1.5678798) q[1];
sx q[1];
rz(-1.5710305) q[1];
x q[2];
rz(-0.51497634) q[3];
sx q[3];
rz(-1.4287717) q[3];
sx q[3];
rz(-1.2339301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2703209) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(-0.16066571) q[2];
rz(1.2284944) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-2.9401949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(-1.6231712) q[1];
sx q[1];
rz(-2.7802614) q[1];
sx q[1];
rz(-2.8526715) q[1];
rz(3.1186947) q[2];
sx q[2];
rz(-1.8635293) q[2];
sx q[2];
rz(0.40018774) q[2];
rz(-2.992441) q[3];
sx q[3];
rz(-2.0059398) q[3];
sx q[3];
rz(-2.5930349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
