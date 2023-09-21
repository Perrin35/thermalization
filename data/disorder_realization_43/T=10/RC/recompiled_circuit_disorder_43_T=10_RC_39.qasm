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
rz(1.8068846) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87603509) q[0];
sx q[0];
rz(-1.8119438) q[0];
sx q[0];
rz(0.59224706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1337778) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(2.9993338) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.098903499) q[1];
sx q[1];
rz(-2.1747723) q[1];
sx q[1];
rz(-0.16278111) q[1];
x q[2];
rz(-2.3372041) q[3];
sx q[3];
rz(-1.8612923) q[3];
sx q[3];
rz(-2.1936072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.3403085) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(0.006342412) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(2.7768163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22390631) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(2.9427337) q[1];
x q[2];
rz(-1.3199602) q[3];
sx q[3];
rz(-2.1187966) q[3];
sx q[3];
rz(2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(-3.045851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11287963) q[0];
sx q[0];
rz(-0.23447795) q[0];
sx q[0];
rz(-2.2857091) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0289621) q[2];
sx q[2];
rz(-2.0320227) q[2];
sx q[2];
rz(1.065965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0248191) q[1];
sx q[1];
rz(-1.0305335) q[1];
sx q[1];
rz(-2.5953672) q[1];
rz(-0.5511958) q[3];
sx q[3];
rz(-1.5407729) q[3];
sx q[3];
rz(2.4195645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948807) q[0];
sx q[0];
rz(-2.046642) q[0];
sx q[0];
rz(-2.4701719) q[0];
rz(-2.1448574) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(0.10276375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32767195) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(3.0753067) q[1];
rz(1.8978118) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(3.0781086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610883) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(-0.9719406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71530576) q[2];
sx q[2];
rz(-1.5537405) q[2];
sx q[2];
rz(-2.7280083) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.638006) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(1.2739146) q[1];
x q[2];
rz(-2.2361034) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(-2.820435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(0.73227698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2039316) q[0];
sx q[0];
rz(-1.7294149) q[0];
sx q[0];
rz(-0.21477867) q[0];
x q[1];
rz(-0.76262577) q[2];
sx q[2];
rz(-2.0848993) q[2];
sx q[2];
rz(-2.5800173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1440891) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(1.287582) q[1];
x q[2];
rz(-0.5703339) q[3];
sx q[3];
rz(-1.3864281) q[3];
sx q[3];
rz(1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(-2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(0.77004534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8806481) q[0];
sx q[0];
rz(-0.62793193) q[0];
sx q[0];
rz(-0.29675608) q[0];
rz(0.99880481) q[2];
sx q[2];
rz(-2.4704128) q[2];
sx q[2];
rz(-0.83516781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3502096) q[1];
sx q[1];
rz(-1.0648799) q[1];
sx q[1];
rz(-0.70178589) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9023864) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(-0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-2.512152) q[2];
sx q[2];
rz(-1.9295613) q[2];
sx q[2];
rz(1.7164873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(-1.3543345) q[1];
x q[2];
rz(3.0931926) q[3];
sx q[3];
rz(-1.0194922) q[3];
sx q[3];
rz(-1.6325523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(3.1270694) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(0.6706388) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(2.3502137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8732849) q[0];
sx q[0];
rz(-2.377016) q[0];
sx q[0];
rz(2.0386049) q[0];
rz(-pi) q[1];
rz(0.4653761) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(-2.344775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(-2.6507069) q[1];
x q[2];
rz(1.295624) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(0.59797257) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35683435) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(0.16108315) q[0];
x q[1];
rz(1.5153377) q[2];
sx q[2];
rz(-2.7055253) q[2];
sx q[2];
rz(-3.1258294) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8163562) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(-2.7951294) q[1];
x q[2];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.0040141) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(2.3168646) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(0.23562283) q[3];
sx q[3];
rz(-1.0141254) q[3];
sx q[3];
rz(0.49401382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
