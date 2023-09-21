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
rz(4.3918443) q[0];
sx q[0];
rz(10.759486) q[0];
rz(2.788738) q[1];
sx q[1];
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87603509) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(0.59224706) q[0];
rz(-pi) q[1];
rz(1.0693597) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(1.2944348) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9589899) q[1];
sx q[1];
rz(-0.62287736) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-2.3372041) q[3];
sx q[3];
rz(-1.2803004) q[3];
sx q[3];
rz(2.1936072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(0.006342412) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20992499) q[0];
sx q[0];
rz(-1.7010265) q[0];
sx q[0];
rz(-1.3757214) q[0];
x q[1];
rz(-0.52829929) q[2];
sx q[2];
rz(-1.0353147) q[2];
sx q[2];
rz(-2.2250125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.983922) q[1];
sx q[1];
rz(-1.6323286) q[1];
sx q[1];
rz(2.8308949) q[1];
rz(-2.7553495) q[3];
sx q[3];
rz(-0.59730232) q[3];
sx q[3];
rz(-2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94770849) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(-0.01097824) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-3.045851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75671065) q[0];
sx q[0];
rz(-1.4178935) q[0];
sx q[0];
rz(-1.3923313) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.1635457) q[2];
sx q[2];
rz(-2.4207052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15624554) q[1];
sx q[1];
rz(-2.3932082) q[1];
sx q[1];
rz(0.85703874) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.057283244) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-0.37477469) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-2.1077164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(0.98866776) q[0];
x q[1];
rz(2.6150377) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(-1.9499792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2305206) q[1];
sx q[1];
rz(-1.5057179) q[1];
sx q[1];
rz(-1.7622403) q[1];
x q[2];
rz(-1.8978118) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(-3.0781086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60526472) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(-1.3461793) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1155903) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(1.1375839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50358665) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(1.867678) q[1];
rz(-pi) q[2];
rz(1.8168713) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(1.0567997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-0.73227698) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4740144) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(-1.7330806) q[0];
rz(-0.68533021) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(2.6076917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5413943) q[1];
sx q[1];
rz(-1.3458034) q[1];
sx q[1];
rz(-0.66585983) q[1];
x q[2];
rz(2.5712588) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(-1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5724065) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(-0.77004534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89966398) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(1.7799737) q[0];
rz(0.40600834) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(-1.6183491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5321977) q[1];
sx q[1];
rz(-0.97071338) q[1];
sx q[1];
rz(-0.9432015) q[1];
rz(-pi) q[2];
rz(0.50290147) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(0.46616947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-0.24169895) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241656) q[0];
sx q[0];
rz(-0.96091849) q[0];
sx q[0];
rz(-0.92745552) q[0];
x q[1];
rz(2.5744152) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(0.30356193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.061325039) q[1];
sx q[1];
rz(-2.419201) q[1];
sx q[1];
rz(-0.25218833) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(-0.036389694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(2.4709539) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489483) q[0];
sx q[0];
rz(-1.2533422) q[0];
sx q[0];
rz(-2.2788458) q[0];
x q[1];
rz(2.6762166) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(-0.79681764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8225122) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(-0.49088571) q[1];
x q[2];
rz(-2.9744042) q[3];
sx q[3];
rz(-1.0382004) q[3];
sx q[3];
rz(-2.5933468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0224232) q[0];
sx q[0];
rz(-1.701231) q[0];
sx q[0];
rz(-2.2020257) q[0];
x q[1];
rz(3.1157687) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(-0.045407427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41439357) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(-0.90331932) q[1];
rz(2.4139666) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(2.0000134) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0040141) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(-0.23562283) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
