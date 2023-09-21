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
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87603509) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(-0.59224706) q[0];
rz(2.8085254) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.8884459) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5648956) q[1];
sx q[1];
rz(-1.4370059) q[1];
sx q[1];
rz(2.1810075) q[1];
rz(0.8043886) q[3];
sx q[3];
rz(-1.8612923) q[3];
sx q[3];
rz(-2.1936072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-0.78380084) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.8012841) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77942383) q[0];
sx q[0];
rz(-2.9075025) q[0];
sx q[0];
rz(0.97658821) q[0];
rz(-2.1727174) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.7768163) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(-2.9427337) q[1];
x q[2];
rz(-0.38624318) q[3];
sx q[3];
rz(-0.59730232) q[3];
sx q[3];
rz(2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851615) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(-0.01097824) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(3.045851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.384882) q[0];
sx q[0];
rz(-1.4178935) q[0];
sx q[0];
rz(1.3923313) q[0];
rz(-2.0289621) q[2];
sx q[2];
rz(-1.1095699) q[2];
sx q[2];
rz(-1.065965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2925551) q[1];
sx q[1];
rz(-2.0325066) q[1];
sx q[1];
rz(-2.1828116) q[1];
rz(-pi) q[2];
rz(-2.5903969) q[3];
sx q[3];
rz(-1.5407729) q[3];
sx q[3];
rz(0.72202819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(2.1077164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046712) q[0];
sx q[0];
rz(-1.0949507) q[0];
sx q[0];
rz(-2.4701719) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3029762) q[2];
sx q[2];
rz(-2.4258483) q[2];
sx q[2];
rz(2.0640304) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32767195) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(3.0753067) q[1];
rz(-pi) q[2];
rz(-2.8598966) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(-1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60526472) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(1.3461793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71530576) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(0.41358435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9628323) q[1];
sx q[1];
rz(-1.8492336) q[1];
sx q[1];
rz(2.7774485) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90548924) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(-2.820435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(2.4093157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4740144) q[0];
sx q[0];
rz(-1.3587553) q[0];
sx q[0];
rz(-1.7330806) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2339091) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(1.6931319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4477168) q[1];
sx q[1];
rz(-0.69732053) q[1];
sx q[1];
rz(0.35481528) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3527649) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(0.41527173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-0.77004534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89966398) q[0];
sx q[0];
rz(-2.1673492) q[0];
sx q[0];
rz(1.7799737) q[0];
rz(-2.1427878) q[2];
sx q[2];
rz(-2.4704128) q[2];
sx q[2];
rz(2.3064248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3502096) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-0.70178589) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(-0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-0.24169895) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0057356) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(2.4321796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.512152) q[2];
sx q[2];
rz(-1.9295613) q[2];
sx q[2];
rz(1.4251054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8230799) q[1];
sx q[1];
rz(-1.4050583) q[1];
sx q[1];
rz(2.4351099) q[1];
rz(2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7991379) q[0];
sx q[0];
rz(-2.236811) q[0];
sx q[0];
rz(0.40823437) q[0];
rz(-pi) q[1];
rz(-0.80957885) q[2];
sx q[2];
rz(-2.5123345) q[2];
sx q[2];
rz(1.4776243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-2.0171595) q[1];
sx q[1];
rz(-1.1091713) q[1];
rz(1.0320372) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62771195) q[0];
sx q[0];
rz(-0.6427592) q[0];
sx q[0];
rz(-1.3520665) q[0];
rz(-0.025823921) q[2];
sx q[2];
rz(-2.0061473) q[2];
sx q[2];
rz(0.045407427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41439357) q[1];
sx q[1];
rz(-2.6135923) q[1];
sx q[1];
rz(-2.2382733) q[1];
rz(-pi) q[2];
rz(-0.72762604) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(2.805368) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
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
rz(-2.3168646) q[2];
sx q[2];
rz(-1.3168954) q[2];
sx q[2];
rz(-1.6864824) q[2];
rz(1.0014793) q[3];
sx q[3];
rz(-1.3712728) q[3];
sx q[3];
rz(-1.2029592) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
