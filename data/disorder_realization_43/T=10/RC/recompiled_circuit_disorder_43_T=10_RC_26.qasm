OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061493) q[0];
sx q[0];
rz(-0.9978928) q[0];
sx q[0];
rz(1.8589622) q[0];
rz(-pi) q[1];
rz(-2.8085254) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.2531467) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5648956) q[1];
sx q[1];
rz(-1.4370059) q[1];
sx q[1];
rz(2.1810075) q[1];
rz(-pi) q[2];
rz(0.39335143) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-0.006342412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(0.13271876) q[0];
rz(2.1727174) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(2.7768163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7087242) q[1];
sx q[1];
rz(-1.2607062) q[1];
sx q[1];
rz(-1.6354145) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56224058) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(-1.7052887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.384882) q[0];
sx q[0];
rz(-1.7236992) q[0];
sx q[0];
rz(1.7492613) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.1635457) q[2];
sx q[2];
rz(-2.4207052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2925551) q[1];
sx q[1];
rz(-2.0325066) q[1];
sx q[1];
rz(-0.95878102) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6060353) q[3];
sx q[3];
rz(-1.0198776) q[3];
sx q[3];
rz(2.2743724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-2.3068008) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9469706) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(2.7950177) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(-2.1529249) q[0];
x q[1];
rz(-0.99673523) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(-0.10276375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8139207) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(-3.0753067) q[1];
rz(-2.8598966) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(1.4178993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(0.57410747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87724553) q[0];
sx q[0];
rz(-2.1640722) q[0];
sx q[0];
rz(-0.1546774) q[0];
rz(-pi) q[1];
rz(-2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(2.7280083) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50358665) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(-1.2739146) q[1];
x q[2];
rz(2.9470109) q[3];
sx q[3];
rz(-0.91472018) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0064272881) q[0];
sx q[0];
rz(-0.26627243) q[0];
sx q[0];
rz(2.4977495) q[0];
x q[1];
rz(-0.90768355) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(1.4484608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69387586) q[1];
sx q[1];
rz(-2.4442721) q[1];
sx q[1];
rz(-0.35481528) q[1];
rz(-pi) q[2];
rz(-2.5712588) q[3];
sx q[3];
rz(-1.3864281) q[3];
sx q[3];
rz(2.1031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89966398) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(-1.7799737) q[0];
rz(0.40600834) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(-1.5232435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5321977) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(-2.1983912) q[1];
rz(-pi) q[2];
rz(-0.55925925) q[3];
sx q[3];
rz(-0.57563215) q[3];
sx q[3];
rz(1.5534793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057356) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(0.70941305) q[0];
x q[1];
rz(-2.0051458) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(0.39603147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(1.7872582) q[1];
rz(-pi) q[2];
rz(1.6493158) q[3];
sx q[3];
rz(-0.55320569) q[3];
sx q[3];
rz(-1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-0.79137897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489483) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(2.2788458) q[0];
rz(-pi) q[1];
rz(2.6762166) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(-0.79681764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.040175166) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(-1.1091713) q[1];
rz(-1.295624) q[3];
sx q[3];
rz(-0.55579805) q[3];
sx q[3];
rz(-0.2273699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(1.3806608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.3520665) q[0];
rz(1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.5362816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41439357) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(-0.90331932) q[1];
rz(-pi) q[2];
rz(-2.9525083) q[3];
sx q[3];
rz(-1.7366647) q[3];
sx q[3];
rz(1.7961111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(2.9059698) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];