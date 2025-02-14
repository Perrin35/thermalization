OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5767515) q[0];
sx q[0];
rz(-0.54902005) q[0];
sx q[0];
rz(0.18000552) q[0];
rz(0.89927468) q[1];
sx q[1];
rz(-0.86714309) q[1];
sx q[1];
rz(1.4049621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1464935) q[0];
sx q[0];
rz(-1.6990663) q[0];
sx q[0];
rz(1.487182) q[0];
x q[1];
rz(-1.399171) q[2];
sx q[2];
rz(-1.020198) q[2];
sx q[2];
rz(-2.0623178) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22706341) q[1];
sx q[1];
rz(-1.8957912) q[1];
sx q[1];
rz(1.090828) q[1];
rz(-pi) q[2];
rz(-1.8990535) q[3];
sx q[3];
rz(-2.0363288) q[3];
sx q[3];
rz(-1.5256603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0286502) q[2];
sx q[2];
rz(-2.3990227) q[2];
sx q[2];
rz(0.58153018) q[2];
rz(0.90315008) q[3];
sx q[3];
rz(-1.4862458) q[3];
sx q[3];
rz(-0.36710292) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0105932) q[0];
sx q[0];
rz(-1.4073263) q[0];
sx q[0];
rz(2.0043688) q[0];
rz(-1.3251023) q[1];
sx q[1];
rz(-2.4749327) q[1];
sx q[1];
rz(2.4773662) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2720075) q[0];
sx q[0];
rz(-1.7427326) q[0];
sx q[0];
rz(0.26845064) q[0];
x q[1];
rz(-1.4777123) q[2];
sx q[2];
rz(-2.4126894) q[2];
sx q[2];
rz(0.62076694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8926516) q[1];
sx q[1];
rz(-2.5499857) q[1];
sx q[1];
rz(-2.4581562) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8898403) q[3];
sx q[3];
rz(-0.24067146) q[3];
sx q[3];
rz(1.9001719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60435158) q[2];
sx q[2];
rz(-1.1553355) q[2];
sx q[2];
rz(1.381116) q[2];
rz(-0.85754496) q[3];
sx q[3];
rz(-0.92567912) q[3];
sx q[3];
rz(3.082357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.183476) q[0];
sx q[0];
rz(-0.97267946) q[0];
sx q[0];
rz(-0.78829515) q[0];
rz(-2.6745785) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(-0.83795396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4350906) q[0];
sx q[0];
rz(-1.52924) q[0];
sx q[0];
rz(2.6986319) q[0];
rz(-1.1594095) q[2];
sx q[2];
rz(-1.5452048) q[2];
sx q[2];
rz(2.7513663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0744255) q[1];
sx q[1];
rz(-0.29109368) q[1];
sx q[1];
rz(1.3240678) q[1];
x q[2];
rz(0.31354745) q[3];
sx q[3];
rz(-1.6308074) q[3];
sx q[3];
rz(2.0504232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9517842) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(-1.5365441) q[2];
rz(-2.8168081) q[3];
sx q[3];
rz(-1.4835417) q[3];
sx q[3];
rz(1.8713846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.558641) q[0];
sx q[0];
rz(-1.0277717) q[0];
sx q[0];
rz(0.047274832) q[0];
rz(-1.6167697) q[1];
sx q[1];
rz(-2.6998417) q[1];
sx q[1];
rz(-1.5025274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29848443) q[0];
sx q[0];
rz(-1.8445948) q[0];
sx q[0];
rz(-0.47852935) q[0];
rz(-pi) q[1];
rz(2.7418361) q[2];
sx q[2];
rz(-2.6538678) q[2];
sx q[2];
rz(-2.9903811) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2983809) q[1];
sx q[1];
rz(-1.3271558) q[1];
sx q[1];
rz(1.3158821) q[1];
x q[2];
rz(-2.8144275) q[3];
sx q[3];
rz(-1.7594271) q[3];
sx q[3];
rz(1.908055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9105685) q[2];
sx q[2];
rz(-0.39414057) q[2];
sx q[2];
rz(-2.2310889) q[2];
rz(-2.2855811) q[3];
sx q[3];
rz(-1.4166219) q[3];
sx q[3];
rz(3.0205309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1055792) q[0];
sx q[0];
rz(-1.7511022) q[0];
sx q[0];
rz(-0.068583071) q[0];
rz(-2.4225281) q[1];
sx q[1];
rz(-1.9870575) q[1];
sx q[1];
rz(2.7209435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7375355) q[0];
sx q[0];
rz(-2.1348663) q[0];
sx q[0];
rz(-0.89667908) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22847036) q[2];
sx q[2];
rz(-1.4632311) q[2];
sx q[2];
rz(3.1247849) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44608403) q[1];
sx q[1];
rz(-0.38818103) q[1];
sx q[1];
rz(-3.0384427) q[1];
rz(0.64780541) q[3];
sx q[3];
rz(-1.2820377) q[3];
sx q[3];
rz(1.6132068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61578304) q[2];
sx q[2];
rz(-2.1685648) q[2];
sx q[2];
rz(-0.77965492) q[2];
rz(-3.0443794) q[3];
sx q[3];
rz(-1.3262871) q[3];
sx q[3];
rz(1.1965082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.6355316) q[0];
sx q[0];
rz(-0.66672915) q[0];
sx q[0];
rz(0.42688236) q[0];
rz(1.6905258) q[1];
sx q[1];
rz(-1.1208813) q[1];
sx q[1];
rz(2.5371187) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647897) q[0];
sx q[0];
rz(-1.2493629) q[0];
sx q[0];
rz(3.1040756) q[0];
rz(0.25570095) q[2];
sx q[2];
rz(-1.2947645) q[2];
sx q[2];
rz(2.7118341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0591914) q[1];
sx q[1];
rz(-1.7520604) q[1];
sx q[1];
rz(3.0742765) q[1];
rz(-pi) q[2];
rz(1.7981619) q[3];
sx q[3];
rz(-1.0195838) q[3];
sx q[3];
rz(2.3662596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84588593) q[2];
sx q[2];
rz(-0.87575951) q[2];
sx q[2];
rz(-2.8952944) q[2];
rz(2.1379499) q[3];
sx q[3];
rz(-1.4469888) q[3];
sx q[3];
rz(0.083855696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698782) q[0];
sx q[0];
rz(-3.0666252) q[0];
sx q[0];
rz(2.9366034) q[0];
rz(0.21944731) q[1];
sx q[1];
rz(-1.4402025) q[1];
sx q[1];
rz(-2.8622389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49528564) q[0];
sx q[0];
rz(-1.4581465) q[0];
sx q[0];
rz(1.1239284) q[0];
rz(-pi) q[1];
rz(1.749239) q[2];
sx q[2];
rz(-2.8861336) q[2];
sx q[2];
rz(-2.5551772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6705888) q[1];
sx q[1];
rz(-1.6274954) q[1];
sx q[1];
rz(0.55135552) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8362813) q[3];
sx q[3];
rz(-2.0274912) q[3];
sx q[3];
rz(1.8165464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67639914) q[2];
sx q[2];
rz(-0.78130829) q[2];
sx q[2];
rz(-2.264273) q[2];
rz(1.8026132) q[3];
sx q[3];
rz(-0.78290144) q[3];
sx q[3];
rz(1.3907998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.029123) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(2.8068338) q[0];
rz(-2.8107457) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(-0.61029339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7772737) q[0];
sx q[0];
rz(-0.36812447) q[0];
sx q[0];
rz(2.0859531) q[0];
rz(2.2525983) q[2];
sx q[2];
rz(-0.51883139) q[2];
sx q[2];
rz(0.31457043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6107294) q[1];
sx q[1];
rz(-1.7708774) q[1];
sx q[1];
rz(-2.6310573) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1842347) q[3];
sx q[3];
rz(-1.6803553) q[3];
sx q[3];
rz(0.027519634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94157964) q[2];
sx q[2];
rz(-2.104685) q[2];
sx q[2];
rz(1.3882136) q[2];
rz(0.16418223) q[3];
sx q[3];
rz(-1.0378342) q[3];
sx q[3];
rz(0.29844704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.91069094) q[0];
sx q[0];
rz(-2.0257484) q[0];
sx q[0];
rz(-0.13634613) q[0];
rz(0.048010437) q[1];
sx q[1];
rz(-2.127357) q[1];
sx q[1];
rz(1.8528574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906517) q[0];
sx q[0];
rz(-1.3539292) q[0];
sx q[0];
rz(1.2614388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1196278) q[2];
sx q[2];
rz(-1.3578312) q[2];
sx q[2];
rz(-2.1599744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9933208) q[1];
sx q[1];
rz(-2.848756) q[1];
sx q[1];
rz(-0.50334986) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39548042) q[3];
sx q[3];
rz(-2.2321475) q[3];
sx q[3];
rz(1.9868074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5137704) q[2];
sx q[2];
rz(-0.0063535293) q[2];
sx q[2];
rz(0.18056907) q[2];
rz(2.0020961) q[3];
sx q[3];
rz(-1.6264911) q[3];
sx q[3];
rz(0.12038055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9554431) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(1.0261616) q[0];
rz(1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(-0.06591448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.41734) q[0];
sx q[0];
rz(-0.98729649) q[0];
sx q[0];
rz(-2.3704607) q[0];
rz(-pi) q[1];
rz(-1.3990551) q[2];
sx q[2];
rz(-1.2615924) q[2];
sx q[2];
rz(-1.1698674) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38635284) q[1];
sx q[1];
rz(-0.59659472) q[1];
sx q[1];
rz(-0.42930023) q[1];
rz(-pi) q[2];
rz(-0.58760175) q[3];
sx q[3];
rz(-1.3580702) q[3];
sx q[3];
rz(-1.0472681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5790448) q[2];
sx q[2];
rz(-1.0479835) q[2];
sx q[2];
rz(2.6143796) q[2];
rz(-0.15050091) q[3];
sx q[3];
rz(-0.42948693) q[3];
sx q[3];
rz(-2.785717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81446205) q[0];
sx q[0];
rz(-0.7939864) q[0];
sx q[0];
rz(-3.0891147) q[0];
rz(1.7809226) q[1];
sx q[1];
rz(-0.7285898) q[1];
sx q[1];
rz(1.8026112) q[1];
rz(0.32329516) q[2];
sx q[2];
rz(-1.9920762) q[2];
sx q[2];
rz(2.6270234) q[2];
rz(0.94447847) q[3];
sx q[3];
rz(-1.2068268) q[3];
sx q[3];
rz(-2.5218388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
