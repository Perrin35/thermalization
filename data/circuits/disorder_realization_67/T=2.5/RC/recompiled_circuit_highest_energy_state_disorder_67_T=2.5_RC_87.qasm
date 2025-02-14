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
rz(1.2713852) q[0];
sx q[0];
rz(-0.013590824) q[0];
sx q[0];
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(-1.3807715) q[1];
sx q[1];
rz(0.071579054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893104) q[0];
sx q[0];
rz(-1.5610862) q[0];
sx q[0];
rz(-0.025061314) q[0];
rz(-pi) q[1];
rz(-2.6144892) q[2];
sx q[2];
rz(-2.6968535) q[2];
sx q[2];
rz(1.2232194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4178582) q[1];
sx q[1];
rz(-1.589554) q[1];
sx q[1];
rz(0.0049519227) q[1];
x q[2];
rz(0.79238331) q[3];
sx q[3];
rz(-0.43607831) q[3];
sx q[3];
rz(2.53873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2804395) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(0.46671483) q[2];
rz(-2.6672582) q[3];
sx q[3];
rz(-0.021183906) q[3];
sx q[3];
rz(-3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83653432) q[0];
sx q[0];
rz(-2.6439522) q[0];
sx q[0];
rz(0.019813892) q[0];
rz(-1.5927947) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(-1.6682495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9112355) q[0];
sx q[0];
rz(-0.94382554) q[0];
sx q[0];
rz(-2.4751653) q[0];
rz(-pi) q[1];
rz(0.25924087) q[2];
sx q[2];
rz(-1.2851641) q[2];
sx q[2];
rz(0.91773568) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2624609) q[1];
sx q[1];
rz(-1.3768798) q[1];
sx q[1];
rz(0.078207774) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4056689) q[3];
sx q[3];
rz(-2.1659327) q[3];
sx q[3];
rz(-1.6079966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8389429) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(1.2897276) q[2];
rz(1.9047811) q[3];
sx q[3];
rz(-2.808282) q[3];
sx q[3];
rz(-0.70241565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692093) q[0];
sx q[0];
rz(-1.1595668) q[0];
sx q[0];
rz(-1.5709391) q[0];
rz(-1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(-2.6872046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4226625) q[0];
sx q[0];
rz(-0.61394982) q[0];
sx q[0];
rz(-1.802868) q[0];
x q[1];
rz(3.1151388) q[2];
sx q[2];
rz(-1.4288386) q[2];
sx q[2];
rz(-0.96645025) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2889298) q[1];
sx q[1];
rz(-1.689012) q[1];
sx q[1];
rz(3.1187936) q[1];
x q[2];
rz(2.347499) q[3];
sx q[3];
rz(-2.3183868) q[3];
sx q[3];
rz(2.2633207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4001974) q[2];
sx q[2];
rz(-3.1252842) q[2];
sx q[2];
rz(-2.7941217) q[2];
rz(-0.45626429) q[3];
sx q[3];
rz(-3.1268692) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149813) q[0];
sx q[0];
rz(-1.240629) q[0];
sx q[0];
rz(1.4062784) q[0];
rz(0.44334385) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(1.571507) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31440526) q[0];
sx q[0];
rz(-1.6256204) q[0];
sx q[0];
rz(2.501295) q[0];
x q[1];
rz(-2.3894044) q[2];
sx q[2];
rz(-0.094399422) q[2];
sx q[2];
rz(2.029325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6199377) q[1];
sx q[1];
rz(-1.8490377) q[1];
sx q[1];
rz(-3.131429) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27914234) q[3];
sx q[3];
rz(-1.1542392) q[3];
sx q[3];
rz(-3.109351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7683679) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(-0.079252871) q[2];
rz(1.9521889) q[3];
sx q[3];
rz(-1.3711843) q[3];
sx q[3];
rz(-1.5090548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.70334148) q[0];
sx q[0];
rz(-2.6282613) q[0];
sx q[0];
rz(-0.80605036) q[0];
rz(2.3015859) q[1];
sx q[1];
rz(-3.1286616) q[1];
sx q[1];
rz(2.3654225) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7568593) q[0];
sx q[0];
rz(-2.1549468) q[0];
sx q[0];
rz(1.7279051) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1311536) q[2];
sx q[2];
rz(-1.569328) q[2];
sx q[2];
rz(-1.8515585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35368109) q[1];
sx q[1];
rz(-3.0056142) q[1];
sx q[1];
rz(-0.087894364) q[1];
x q[2];
rz(-0.40078409) q[3];
sx q[3];
rz(-2.8607607) q[3];
sx q[3];
rz(0.1268445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.025803056) q[2];
sx q[2];
rz(-1.5610521) q[2];
sx q[2];
rz(-0.72581327) q[2];
rz(0.18802655) q[3];
sx q[3];
rz(-0.060421061) q[3];
sx q[3];
rz(-2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96227729) q[0];
sx q[0];
rz(-2.582452) q[0];
sx q[0];
rz(0.54139262) q[0];
rz(-2.9560282) q[1];
sx q[1];
rz(-1.5488397) q[1];
sx q[1];
rz(-0.12841368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1379515) q[0];
sx q[0];
rz(-1.9310474) q[0];
sx q[0];
rz(-1.2482578) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1368581) q[2];
sx q[2];
rz(-1.4494697) q[2];
sx q[2];
rz(1.5743299) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0182788) q[1];
sx q[1];
rz(-1.6524846) q[1];
sx q[1];
rz(-0.42079349) q[1];
rz(-pi) q[2];
rz(1.382368) q[3];
sx q[3];
rz(-1.5811833) q[3];
sx q[3];
rz(-2.5760108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7554756) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(-0.84121394) q[2];
rz(-0.20127131) q[3];
sx q[3];
rz(-1.5320675) q[3];
sx q[3];
rz(-2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.5801308) q[0];
sx q[0];
rz(-2.353297) q[0];
sx q[0];
rz(-1.5607675) q[0];
rz(-0.67281094) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(3.1127081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88850856) q[0];
sx q[0];
rz(-2.7612503) q[0];
sx q[0];
rz(2.4230291) q[0];
rz(-0.47841623) q[2];
sx q[2];
rz(-1.3923651) q[2];
sx q[2];
rz(-0.24716864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4116482) q[1];
sx q[1];
rz(-1.8081417) q[1];
sx q[1];
rz(1.9945595) q[1];
x q[2];
rz(-0.77856346) q[3];
sx q[3];
rz(-1.1751919) q[3];
sx q[3];
rz(2.9958519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9081356) q[2];
sx q[2];
rz(-2.5448749) q[2];
sx q[2];
rz(-0.92823589) q[2];
rz(-2.8440031) q[3];
sx q[3];
rz(-0.15407763) q[3];
sx q[3];
rz(-2.2969864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0723648) q[0];
sx q[0];
rz(-2.897825) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(-2.15436) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(-0.5388906) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19463704) q[0];
sx q[0];
rz(-1.4879003) q[0];
sx q[0];
rz(3.1350053) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1168836) q[2];
sx q[2];
rz(-1.5214143) q[2];
sx q[2];
rz(2.1076237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1325593) q[1];
sx q[1];
rz(-1.0485639) q[1];
sx q[1];
rz(-1.9886677) q[1];
x q[2];
rz(1.5032926) q[3];
sx q[3];
rz(-1.5933523) q[3];
sx q[3];
rz(-1.8543138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99523669) q[2];
sx q[2];
rz(-3.1260999) q[2];
sx q[2];
rz(2.8323979) q[2];
rz(2.501798) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(3.0777847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(0.30534202) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(-3.0869361) q[0];
rz(-2.0089741) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(1.2861015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9817106) q[0];
sx q[0];
rz(-1.3927841) q[0];
sx q[0];
rz(-0.044248754) q[0];
x q[1];
rz(3.1007721) q[2];
sx q[2];
rz(-1.5287182) q[2];
sx q[2];
rz(0.034860858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1680345) q[1];
sx q[1];
rz(-1.66093) q[1];
sx q[1];
rz(0.22613495) q[1];
x q[2];
rz(2.0580106) q[3];
sx q[3];
rz(-1.6346674) q[3];
sx q[3];
rz(2.4364249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.22141156) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(-2.0378225) q[2];
rz(-1.6086027) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(-2.485763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00103818) q[0];
sx q[0];
rz(-2.9619205) q[0];
sx q[0];
rz(-3.1371064) q[0];
rz(1.5890315) q[1];
sx q[1];
rz(-1.4493425) q[1];
sx q[1];
rz(3.0844614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308459) q[0];
sx q[0];
rz(-1.391469) q[0];
sx q[0];
rz(-0.097436949) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22259076) q[2];
sx q[2];
rz(-1.4302956) q[2];
sx q[2];
rz(2.8369571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88567121) q[1];
sx q[1];
rz(-2.2388865) q[1];
sx q[1];
rz(0.1541962) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7369676) q[3];
sx q[3];
rz(-0.7471841) q[3];
sx q[3];
rz(2.2053444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.750018) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(0.86891437) q[2];
rz(0.9817552) q[3];
sx q[3];
rz(-0.029564094) q[3];
sx q[3];
rz(0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960196) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(-2.6847196) q[1];
sx q[1];
rz(-2.9869106) q[1];
sx q[1];
rz(3.0966495) q[1];
rz(-2.9484684) q[2];
sx q[2];
rz(-0.77342351) q[2];
sx q[2];
rz(0.20062994) q[2];
rz(-0.66309912) q[3];
sx q[3];
rz(-1.6840006) q[3];
sx q[3];
rz(-1.4406289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
