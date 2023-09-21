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
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1057518) q[0];
sx q[0];
rz(-0.63397206) q[0];
sx q[0];
rz(0.41497725) q[0];
rz(-pi) q[1];
rz(2.8085254) q[2];
sx q[2];
rz(-1.0330079) q[2];
sx q[2];
rz(-1.8884459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9589899) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.8014924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39335143) q[3];
sx q[3];
rz(-0.8439807) q[3];
sx q[3];
rz(0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(0.006342412) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(3.0088739) q[0];
x q[1];
rz(-0.86654051) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(-1.7689592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43286846) q[1];
sx q[1];
rz(-1.2607062) q[1];
sx q[1];
rz(1.6354145) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5793521) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-0.095741622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.384882) q[0];
sx q[0];
rz(-1.7236992) q[0];
sx q[0];
rz(-1.7492613) q[0];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(2.4207052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84903753) q[1];
sx q[1];
rz(-2.0325066) q[1];
sx q[1];
rz(-2.1828116) q[1];
rz(0.057283244) q[3];
sx q[3];
rz(-0.55192845) q[3];
sx q[3];
rz(2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0694971) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(2.7950177) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046712) q[0];
sx q[0];
rz(-2.046642) q[0];
sx q[0];
rz(-0.67142077) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83861645) q[2];
sx q[2];
rz(-2.4258483) q[2];
sx q[2];
rz(2.0640304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66400601) q[1];
sx q[1];
rz(-0.2020745) q[1];
sx q[1];
rz(1.2408153) q[1];
rz(-2.8598966) q[3];
sx q[3];
rz(-1.2558189) q[3];
sx q[3];
rz(-1.4178993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610883) q[0];
sx q[0];
rz(-1.6988806) q[0];
sx q[0];
rz(-2.1696521) q[0];
x q[1];
rz(2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(-2.7280083) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50358665) q[1];
sx q[1];
rz(-1.9203016) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3247213) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(-1.0567997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(-2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-0.73227698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(1.408512) q[0];
rz(0.90768355) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(1.6931319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1440891) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-pi) q[2];
rz(-0.332571) q[3];
sx q[3];
rz(-0.59623527) q[3];
sx q[3];
rz(-0.81070825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-0.77004534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609445) q[0];
sx q[0];
rz(-0.62793193) q[0];
sx q[0];
rz(0.29675608) q[0];
rz(-0.40600834) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(-1.6183491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7913831) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-2.4398068) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5823334) q[3];
sx q[3];
rz(-0.57563215) q[3];
sx q[3];
rz(-1.5881133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05224932) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(0.71787562) q[0];
x q[1];
rz(-0.56717746) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(0.30356193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0802676) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(2.8894043) q[1];
rz(1.0189692) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(2.3502137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34245472) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(-2.7333583) q[0];
rz(-pi) q[1];
rz(2.0558526) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(0.56318356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(-0.49088571) q[1];
x q[2];
rz(-1.0320372) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(1.1080351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.3806608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1191694) q[0];
sx q[0];
rz(-1.4403617) q[0];
sx q[0];
rz(2.2020257) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1157687) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(0.045407427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8163562) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(0.34646323) q[1];
rz(-2.4139666) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
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
rz(2.1401134) q[3];
sx q[3];
rz(-1.7703198) q[3];
sx q[3];
rz(1.9386335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
