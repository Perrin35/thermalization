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
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655576) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(-0.59224706) q[0];
rz(-2.8085254) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.2531467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.576697) q[1];
sx q[1];
rz(-1.7045867) q[1];
sx q[1];
rz(-0.96058515) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9777771) q[3];
sx q[3];
rz(-0.80899901) q[3];
sx q[3];
rz(0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(0.006342412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7550678) q[0];
sx q[0];
rz(-1.3773943) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(-0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(2.7768163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22390631) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(-2.9427337) q[1];
rz(1.8216324) q[3];
sx q[3];
rz(-2.1187966) q[3];
sx q[3];
rz(2.8743924) q[3];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84155267) q[0];
sx q[0];
rz(-1.3944355) q[0];
sx q[0];
rz(-2.9862613) q[0];
rz(-pi) q[1];
rz(0.7272561) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(-0.22932316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0248191) q[1];
sx q[1];
rz(-1.0305335) q[1];
sx q[1];
rz(-0.54622548) q[1];
rz(3.0843094) q[3];
sx q[3];
rz(-0.55192845) q[3];
sx q[3];
rz(0.79997593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9469706) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82499408) q[0];
sx q[0];
rz(-2.1567417) q[0];
sx q[0];
rz(-2.1529249) q[0];
rz(-0.526555) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(-1.9499792) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66400601) q[1];
sx q[1];
rz(-2.9395182) q[1];
sx q[1];
rz(1.2408153) q[1];
rz(1.2437808) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.7948077) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(0.57410747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60526472) q[0];
sx q[0];
rz(-2.5308373) q[0];
sx q[0];
rz(1.3461793) q[0];
rz(0.026002361) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(1.1375839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9628323) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(2.7774485) q[1];
rz(-pi) q[2];
rz(-2.2361034) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(-2.820435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(-0.69818991) q[0];
rz(-2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-2.4093157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0064272881) q[0];
sx q[0];
rz(-2.8753202) q[0];
sx q[0];
rz(0.64384319) q[0];
rz(-pi) q[1];
rz(-2.3789669) q[2];
sx q[2];
rz(-2.0848993) q[2];
sx q[2];
rz(-0.56157535) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5413943) q[1];
sx q[1];
rz(-1.7957893) q[1];
sx q[1];
rz(-2.4757328) q[1];
x q[2];
rz(-2.8090217) q[3];
sx q[3];
rz(-2.5453574) q[3];
sx q[3];
rz(-0.81070825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-0.15643315) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(0.77004534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8806481) q[0];
sx q[0];
rz(-2.5136607) q[0];
sx q[0];
rz(2.8448366) q[0];
x q[1];
rz(-2.7355843) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(1.6183491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3502096) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(2.4398068) q[1];
rz(1.9023864) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(-2.1945206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(-2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05224932) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(0.71787562) q[0];
rz(-pi) q[1];
rz(2.512152) q[2];
sx q[2];
rz(-1.2120314) q[2];
sx q[2];
rz(1.7164873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.3543345) q[1];
rz(-pi) q[2];
rz(-2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.83207399) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34245472) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(-2.7333583) q[0];
rz(-pi) q[1];
rz(-1.0857401) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(-2.5784091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2453354) q[1];
sx q[1];
rz(-0.63071139) q[1];
sx q[1];
rz(0.74951042) q[1];
rz(-pi) q[2];
rz(-2.1095554) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(2.5436201) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.7609319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1191694) q[0];
sx q[0];
rz(-1.4403617) q[0];
sx q[0];
rz(2.2020257) q[0];
rz(-pi) q[1];
rz(-0.025823921) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(3.0961852) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32523649) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(-2.7951294) q[1];
rz(-1.7396183) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(2.805368) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
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
rz(-1.2119157) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
