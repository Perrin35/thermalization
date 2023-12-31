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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358409) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(2.7266154) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.072233) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(-1.8471579) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5648956) q[1];
sx q[1];
rz(-1.7045867) q[1];
sx q[1];
rz(2.1810075) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7482412) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
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
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(0.86654051) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(1.7689592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.43286846) q[1];
sx q[1];
rz(-1.8808865) q[1];
sx q[1];
rz(-1.5061782) q[1];
x q[2];
rz(-0.38624318) q[3];
sx q[3];
rz(-0.59730232) q[3];
sx q[3];
rz(-0.18988767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-3.1356964) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75671065) q[0];
sx q[0];
rz(-1.4178935) q[0];
sx q[0];
rz(1.3923313) q[0];
rz(1.1126306) q[2];
sx q[2];
rz(-1.1095699) q[2];
sx q[2];
rz(-1.065965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2925551) q[1];
sx q[1];
rz(-1.109086) q[1];
sx q[1];
rz(2.1828116) q[1];
x q[2];
rz(1.6060353) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(-0.86722022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(-0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(0.52571458) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(1.0338763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(0.69181504) q[0];
rz(2.6150377) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(-1.1916135) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66400601) q[1];
sx q[1];
rz(-0.2020745) q[1];
sx q[1];
rz(-1.9007773) q[1];
x q[2];
rz(-1.2437808) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088257) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78050437) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(0.9719406) q[0];
rz(-0.026002361) q[2];
sx q[2];
rz(-0.71547316) q[2];
sx q[2];
rz(1.1375839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9628323) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(0.36414418) q[1];
rz(-1.3247213) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(-2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-0.73227698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(2.926814) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2339091) q[2];
sx q[2];
rz(-2.216202) q[2];
sx q[2];
rz(-1.4484608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4477168) q[1];
sx q[1];
rz(-2.4442721) q[1];
sx q[1];
rz(2.7867774) q[1];
x q[2];
rz(0.5703339) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(-2.1031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-0.15643315) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(-2.3715473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55243385) q[0];
sx q[0];
rz(-1.3981515) q[0];
sx q[0];
rz(0.60683672) q[0];
rz(-2.1427878) q[2];
sx q[2];
rz(-0.6711798) q[2];
sx q[2];
rz(0.83516781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29994446) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(2.4323835) q[1];
x q[2];
rz(-0.55925925) q[3];
sx q[3];
rz(-0.57563215) q[3];
sx q[3];
rz(-1.5881133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05224932) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(2.423717) q[0];
rz(-2.0051458) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(0.39603147) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.061325039) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(0.25218833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0189692) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(0.86274685) q[0];
rz(-pi) q[1];
rz(-0.4653761) q[2];
sx q[2];
rz(-1.1306136) q[2];
sx q[2];
rz(-2.344775) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3190805) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(-2.6507069) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0320372) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(-1.1080351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-0.6427592) q[0];
sx q[0];
rz(-1.7895262) q[0];
rz(-pi) q[1];
rz(-2.0062749) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.5362816) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8163562) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(0.34646323) q[1];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(-2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-2.0000134) q[3];
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
rz(-pi) q[0];
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
rz(0.33967321) q[2];
sx q[2];
rz(-0.85396955) q[2];
sx q[2];
rz(0.11243482) q[2];
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
