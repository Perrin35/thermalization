OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(-1.891341) q[0];
sx q[0];
rz(1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1057518) q[0];
sx q[0];
rz(-0.63397206) q[0];
sx q[0];
rz(-2.7266154) q[0];
rz(-0.33306723) q[2];
sx q[2];
rz(-1.0330079) q[2];
sx q[2];
rz(1.2531467) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18260278) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9777771) q[3];
sx q[3];
rz(-0.80899901) q[3];
sx q[3];
rz(-0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.3403085) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9316677) q[0];
sx q[0];
rz(-1.7010265) q[0];
sx q[0];
rz(-1.7658712) q[0];
x q[1];
rz(-0.52829929) q[2];
sx q[2];
rz(-1.0353147) q[2];
sx q[2];
rz(-2.2250125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(-0.19885893) q[1];
x q[2];
rz(-1.8216324) q[3];
sx q[3];
rz(-1.0227961) q[3];
sx q[3];
rz(2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(-3.1306144) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84155267) q[0];
sx q[0];
rz(-1.7471572) q[0];
sx q[0];
rz(-2.9862613) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6355866) q[2];
sx q[2];
rz(-1.1635457) q[2];
sx q[2];
rz(-2.4207052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15624554) q[1];
sx q[1];
rz(-2.3932082) q[1];
sx q[1];
rz(2.2845539) q[1];
x q[2];
rz(-3.0843094) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(0.79997593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(2.766818) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3165986) q[0];
sx q[0];
rz(-2.1567417) q[0];
sx q[0];
rz(-2.1529249) q[0];
x q[1];
rz(0.526555) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(1.1916135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8139207) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(-3.0753067) q[1];
rz(-pi) q[2];
rz(-2.2771308) q[3];
sx q[3];
rz(-0.41941386) q[3];
sx q[3];
rz(0.97233397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-2.5674852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78050437) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(2.1696521) q[0];
x q[1];
rz(2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(-2.7280083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9628323) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(-0.36414418) q[1];
x q[2];
rz(-2.9470109) q[3];
sx q[3];
rz(-2.2268725) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(-2.926814) q[0];
rz(-pi) q[1];
rz(0.68533021) q[2];
sx q[2];
rz(-0.88980674) q[2];
sx q[2];
rz(2.6076917) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5413943) q[1];
sx q[1];
rz(-1.7957893) q[1];
sx q[1];
rz(2.4757328) q[1];
rz(-pi) q[2];
rz(-0.332571) q[3];
sx q[3];
rz(-2.5453574) q[3];
sx q[3];
rz(0.81070825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2609445) q[0];
sx q[0];
rz(-0.62793193) q[0];
sx q[0];
rz(2.8448366) q[0];
rz(-pi) q[1];
rz(0.98203512) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(-2.8729168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3502096) q[1];
sx q[1];
rz(-1.0648799) q[1];
sx q[1];
rz(2.4398068) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2392063) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1358571) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(2.4321796) q[0];
rz(-0.56717746) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(-2.8380307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3185127) q[1];
sx q[1];
rz(-1.4050583) q[1];
sx q[1];
rz(-0.70648273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.048400064) q[3];
sx q[3];
rz(-1.0194922) q[3];
sx q[3];
rz(-1.6325523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-2.3502137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(0.86274685) q[0];
x q[1];
rz(-2.3320138) q[2];
sx q[2];
rz(-0.62925816) q[2];
sx q[2];
rz(-1.6639683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(2.0324213) q[1];
rz(-2.1095554) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62771195) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(-1.7895262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.626255) q[2];
sx q[2];
rz(-2.7055253) q[2];
sx q[2];
rz(-3.1258294) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8163562) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(-2.7951294) q[1];
rz(-pi) q[2];
rz(-1.7396183) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(-0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(-1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(0.82472807) q[2];
sx q[2];
rz(-1.3168954) q[2];
sx q[2];
rz(-1.6864824) q[2];
rz(-1.929677) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
