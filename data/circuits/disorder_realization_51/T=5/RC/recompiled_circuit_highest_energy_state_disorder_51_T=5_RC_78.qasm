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
rz(0.60740745) q[0];
sx q[0];
rz(-0.4439126) q[0];
sx q[0];
rz(-2.8191415) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(-2.9236887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3716952) q[0];
sx q[0];
rz(-1.2697233) q[0];
sx q[0];
rz(-2.6825344) q[0];
rz(-pi) q[1];
rz(0.18677752) q[2];
sx q[2];
rz(-1.412027) q[2];
sx q[2];
rz(-1.1940317) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17927781) q[1];
sx q[1];
rz(-1.6861818) q[1];
sx q[1];
rz(0.25956599) q[1];
x q[2];
rz(2.3090906) q[3];
sx q[3];
rz(-0.41838405) q[3];
sx q[3];
rz(2.7184021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13088642) q[2];
sx q[2];
rz(-0.43121269) q[2];
sx q[2];
rz(0.19403379) q[2];
rz(-2.7583097) q[3];
sx q[3];
rz(-2.2563939) q[3];
sx q[3];
rz(-1.6398199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5633157) q[0];
sx q[0];
rz(-2.0864154) q[0];
sx q[0];
rz(2.03736) q[0];
rz(-0.68179321) q[1];
sx q[1];
rz(-1.3673404) q[1];
sx q[1];
rz(-1.404748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0007219) q[0];
sx q[0];
rz(-0.28647505) q[0];
sx q[0];
rz(2.6764144) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4787729) q[2];
sx q[2];
rz(-2.347555) q[2];
sx q[2];
rz(-2.8341849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7182285) q[1];
sx q[1];
rz(-1.8830788) q[1];
sx q[1];
rz(0.52257089) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25481635) q[3];
sx q[3];
rz(-1.4597893) q[3];
sx q[3];
rz(-2.1141583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3746987) q[2];
sx q[2];
rz(-2.1465492) q[2];
sx q[2];
rz(-1.9083171) q[2];
rz(-2.342566) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(0.12921216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1610573) q[0];
sx q[0];
rz(-1.072847) q[0];
sx q[0];
rz(-0.55915731) q[0];
rz(-1.1178389) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(-3.0303755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1532899) q[0];
sx q[0];
rz(-0.39228253) q[0];
sx q[0];
rz(-2.8610703) q[0];
rz(-1.8410013) q[2];
sx q[2];
rz(-0.96365721) q[2];
sx q[2];
rz(0.48838136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8168252) q[1];
sx q[1];
rz(-2.2681464) q[1];
sx q[1];
rz(-2.8945385) q[1];
x q[2];
rz(-2.7484557) q[3];
sx q[3];
rz(-1.434978) q[3];
sx q[3];
rz(-2.6640662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8289651) q[2];
sx q[2];
rz(-1.6964922) q[2];
sx q[2];
rz(-1.1980537) q[2];
rz(-1.7184006) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(2.6660582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6296185) q[0];
sx q[0];
rz(-0.29933512) q[0];
sx q[0];
rz(2.3053115) q[0];
rz(1.5760999) q[1];
sx q[1];
rz(-2.0471768) q[1];
sx q[1];
rz(-0.85087585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880575) q[0];
sx q[0];
rz(-1.5374105) q[0];
sx q[0];
rz(-0.99793156) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81340547) q[2];
sx q[2];
rz(-2.4508173) q[2];
sx q[2];
rz(-1.9936313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16164348) q[1];
sx q[1];
rz(-1.4323712) q[1];
sx q[1];
rz(0.78236492) q[1];
rz(-pi) q[2];
rz(0.20765813) q[3];
sx q[3];
rz(-1.1829783) q[3];
sx q[3];
rz(-0.70630276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0389082) q[2];
sx q[2];
rz(-2.2816198) q[2];
sx q[2];
rz(2.6343708) q[2];
rz(2.475259) q[3];
sx q[3];
rz(-2.3907876) q[3];
sx q[3];
rz(0.91608086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0136593) q[0];
sx q[0];
rz(-1.4480696) q[0];
sx q[0];
rz(1.6130945) q[0];
rz(-1.5799892) q[1];
sx q[1];
rz(-1.0630307) q[1];
sx q[1];
rz(-0.16435057) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2737758) q[0];
sx q[0];
rz(-1.651433) q[0];
sx q[0];
rz(0.095726526) q[0];
rz(-pi) q[1];
rz(-0.63173024) q[2];
sx q[2];
rz(-2.1513878) q[2];
sx q[2];
rz(-0.7484865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4867465) q[1];
sx q[1];
rz(-0.89895144) q[1];
sx q[1];
rz(-0.44507546) q[1];
rz(-pi) q[2];
rz(0.23779121) q[3];
sx q[3];
rz(-1.2850518) q[3];
sx q[3];
rz(-1.3956436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5483755) q[2];
sx q[2];
rz(-0.84725738) q[2];
sx q[2];
rz(2.4753921) q[2];
rz(0.14144746) q[3];
sx q[3];
rz(-1.5902767) q[3];
sx q[3];
rz(-0.9945873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87404609) q[0];
sx q[0];
rz(-2.5135437) q[0];
sx q[0];
rz(-2.6603267) q[0];
rz(2.4080243) q[1];
sx q[1];
rz(-1.2679408) q[1];
sx q[1];
rz(3.0487294) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294612) q[0];
sx q[0];
rz(-0.67352897) q[0];
sx q[0];
rz(-0.81443925) q[0];
x q[1];
rz(-1.5687555) q[2];
sx q[2];
rz(-1.5321595) q[2];
sx q[2];
rz(2.3734951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0561547) q[1];
sx q[1];
rz(-2.3533818) q[1];
sx q[1];
rz(-0.5931535) q[1];
rz(0.95488432) q[3];
sx q[3];
rz(-0.43514565) q[3];
sx q[3];
rz(0.26538119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4197454) q[2];
sx q[2];
rz(-1.5200619) q[2];
sx q[2];
rz(-1.9160371) q[2];
rz(-0.078977481) q[3];
sx q[3];
rz(-1.1342528) q[3];
sx q[3];
rz(-1.8592161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52800286) q[0];
sx q[0];
rz(-2.4230175) q[0];
sx q[0];
rz(0.80793107) q[0];
rz(1.5241874) q[1];
sx q[1];
rz(-0.15847358) q[1];
sx q[1];
rz(2.4375516) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2815783) q[0];
sx q[0];
rz(-1.5905955) q[0];
sx q[0];
rz(0.0088282882) q[0];
rz(2.03116) q[2];
sx q[2];
rz(-3.1322479) q[2];
sx q[2];
rz(-1.2197987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0800533) q[1];
sx q[1];
rz(-2.1916917) q[1];
sx q[1];
rz(-0.027959221) q[1];
rz(-pi) q[2];
rz(-2.1857587) q[3];
sx q[3];
rz(-1.4356239) q[3];
sx q[3];
rz(-1.2379902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.070270553) q[2];
sx q[2];
rz(-0.68838781) q[2];
sx q[2];
rz(-2.2570611) q[2];
rz(-0.5136579) q[3];
sx q[3];
rz(-1.9653886) q[3];
sx q[3];
rz(-0.48532143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014548253) q[0];
sx q[0];
rz(-0.7651279) q[0];
sx q[0];
rz(2.2807518) q[0];
rz(-0.035482081) q[1];
sx q[1];
rz(-1.9396962) q[1];
sx q[1];
rz(-1.6532345) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299773) q[0];
sx q[0];
rz(-1.9693144) q[0];
sx q[0];
rz(2.7048955) q[0];
x q[1];
rz(0.22895892) q[2];
sx q[2];
rz(-1.8081852) q[2];
sx q[2];
rz(-1.4942394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1775629) q[1];
sx q[1];
rz(-0.42516758) q[1];
sx q[1];
rz(-1.1336826) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27542563) q[3];
sx q[3];
rz(-1.6622433) q[3];
sx q[3];
rz(2.8851654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52582467) q[2];
sx q[2];
rz(-1.1349698) q[2];
sx q[2];
rz(1.6305249) q[2];
rz(0.75119558) q[3];
sx q[3];
rz(-2.5978751) q[3];
sx q[3];
rz(-0.80120075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649331) q[0];
sx q[0];
rz(-1.7486005) q[0];
sx q[0];
rz(-0.15889731) q[0];
rz(1.502602) q[1];
sx q[1];
rz(-1.331012) q[1];
sx q[1];
rz(1.9749036) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5432271) q[0];
sx q[0];
rz(-2.1413714) q[0];
sx q[0];
rz(-1.8680509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27952898) q[2];
sx q[2];
rz(-1.2567826) q[2];
sx q[2];
rz(1.90998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9433971) q[1];
sx q[1];
rz(-1.5540197) q[1];
sx q[1];
rz(-1.5328963) q[1];
rz(-3.1205584) q[3];
sx q[3];
rz(-1.6004749) q[3];
sx q[3];
rz(-2.6231387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4022973) q[2];
sx q[2];
rz(-0.6604971) q[2];
sx q[2];
rz(2.8759549) q[2];
rz(1.9499251) q[3];
sx q[3];
rz(-1.3375514) q[3];
sx q[3];
rz(-0.55408293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83248508) q[0];
sx q[0];
rz(-2.5256248) q[0];
sx q[0];
rz(-0.68514222) q[0];
rz(2.5849672) q[1];
sx q[1];
rz(-1.4645422) q[1];
sx q[1];
rz(2.305078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292824) q[0];
sx q[0];
rz(-2.5959281) q[0];
sx q[0];
rz(-0.87581234) q[0];
x q[1];
rz(2.0903953) q[2];
sx q[2];
rz(-1.300989) q[2];
sx q[2];
rz(-0.74120159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25364177) q[1];
sx q[1];
rz(-2.3440323) q[1];
sx q[1];
rz(-2.5152999) q[1];
x q[2];
rz(1.8440866) q[3];
sx q[3];
rz(-1.7215842) q[3];
sx q[3];
rz(-2.085872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69540018) q[2];
sx q[2];
rz(-1.7170898) q[2];
sx q[2];
rz(-2.2851473) q[2];
rz(-2.942318) q[3];
sx q[3];
rz(-2.1592185) q[3];
sx q[3];
rz(0.91163951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0675426) q[0];
sx q[0];
rz(-1.0150801) q[0];
sx q[0];
rz(-1.6652921) q[0];
rz(-0.53786565) q[1];
sx q[1];
rz(-1.8712416) q[1];
sx q[1];
rz(1.3457294) q[1];
rz(2.3026148) q[2];
sx q[2];
rz(-0.92409033) q[2];
sx q[2];
rz(1.4983231) q[2];
rz(-2.2945401) q[3];
sx q[3];
rz(-1.3624227) q[3];
sx q[3];
rz(-0.25400454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
