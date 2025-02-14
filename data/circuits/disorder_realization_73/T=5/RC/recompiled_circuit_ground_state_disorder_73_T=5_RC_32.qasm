OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(-2.2247563) q[0];
sx q[0];
rz(0.43489781) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(-2.6225852) q[1];
sx q[1];
rz(2.5595698) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618185) q[0];
sx q[0];
rz(-1.4167182) q[0];
sx q[0];
rz(1.3067354) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2088654) q[2];
sx q[2];
rz(-1.9550394) q[2];
sx q[2];
rz(2.8926433) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4705321) q[1];
sx q[1];
rz(-1.8166421) q[1];
sx q[1];
rz(2.0162986) q[1];
rz(-pi) q[2];
x q[2];
rz(1.01389) q[3];
sx q[3];
rz(-2.8134973) q[3];
sx q[3];
rz(-1.2893639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5883098) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(0.29169875) q[2];
rz(0.45082539) q[3];
sx q[3];
rz(-0.25142938) q[3];
sx q[3];
rz(0.60744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400103) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(1.095358) q[0];
rz(-1.1913242) q[1];
sx q[1];
rz(-0.92000633) q[1];
sx q[1];
rz(-2.2390168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13690925) q[0];
sx q[0];
rz(-2.1666514) q[0];
sx q[0];
rz(-2.3803902) q[0];
rz(-pi) q[1];
rz(0.29400715) q[2];
sx q[2];
rz(-2.2180811) q[2];
sx q[2];
rz(-0.33366007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0863994) q[1];
sx q[1];
rz(-0.31685778) q[1];
sx q[1];
rz(0.39519127) q[1];
rz(-0.91530494) q[3];
sx q[3];
rz(-0.73507959) q[3];
sx q[3];
rz(0.32441586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9537182) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(-0.11889674) q[2];
rz(-0.64905727) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(1.6868748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.651799) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(2.6258262) q[0];
rz(1.0924529) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(1.8663503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37947734) q[0];
sx q[0];
rz(-1.5526875) q[0];
sx q[0];
rz(0.041569592) q[0];
rz(-0.31351201) q[2];
sx q[2];
rz(-1.3517153) q[2];
sx q[2];
rz(-2.8323725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.470388) q[1];
sx q[1];
rz(-0.71031308) q[1];
sx q[1];
rz(1.2673668) q[1];
rz(-0.10026513) q[3];
sx q[3];
rz(-0.79660049) q[3];
sx q[3];
rz(1.0880566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.065993) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(-0.48948151) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(-0.71973962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8664261) q[0];
sx q[0];
rz(-2.8186099) q[0];
sx q[0];
rz(-0.48165709) q[0];
rz(-0.9306759) q[1];
sx q[1];
rz(-1.296867) q[1];
sx q[1];
rz(0.01893386) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1279432) q[0];
sx q[0];
rz(-3.0873723) q[0];
sx q[0];
rz(1.9677866) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90831129) q[2];
sx q[2];
rz(-1.6056955) q[2];
sx q[2];
rz(-2.1878583) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.150731) q[1];
sx q[1];
rz(-0.98349893) q[1];
sx q[1];
rz(2.2022922) q[1];
rz(0.76374028) q[3];
sx q[3];
rz(-2.1472048) q[3];
sx q[3];
rz(-1.9595722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9229752) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(-2.3962928) q[2];
rz(-0.37799147) q[3];
sx q[3];
rz(-1.1443006) q[3];
sx q[3];
rz(-2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.9208263) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(1.1908603) q[0];
rz(0.4885172) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(-0.8078422) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5295083) q[0];
sx q[0];
rz(-2.8017462) q[0];
sx q[0];
rz(-1.8456949) q[0];
rz(2.9132782) q[2];
sx q[2];
rz(-2.0234152) q[2];
sx q[2];
rz(-1.376898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1792308) q[1];
sx q[1];
rz(-0.46474248) q[1];
sx q[1];
rz(-0.19728139) q[1];
rz(-1.8459365) q[3];
sx q[3];
rz(-2.5003331) q[3];
sx q[3];
rz(0.99320557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1390344) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(-0.16415088) q[2];
rz(-0.74603355) q[3];
sx q[3];
rz(-1.371871) q[3];
sx q[3];
rz(3.0677838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1604851) q[0];
sx q[0];
rz(-1.2657413) q[0];
sx q[0];
rz(-2.4142081) q[0];
rz(-2.6630867) q[1];
sx q[1];
rz(-1.452927) q[1];
sx q[1];
rz(2.4047638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664624) q[0];
sx q[0];
rz(-1.7084165) q[0];
sx q[0];
rz(0.055995106) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8520459) q[2];
sx q[2];
rz(-1.6967456) q[2];
sx q[2];
rz(1.151786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0733395) q[1];
sx q[1];
rz(-0.96594772) q[1];
sx q[1];
rz(-2.6920435) q[1];
rz(-pi) q[2];
rz(-3.1276631) q[3];
sx q[3];
rz(-1.1831814) q[3];
sx q[3];
rz(-0.25957169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9937146) q[2];
sx q[2];
rz(-1.0241877) q[2];
sx q[2];
rz(-0.68515879) q[2];
rz(-2.1327175) q[3];
sx q[3];
rz(-2.4515371) q[3];
sx q[3];
rz(0.25079295) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20550263) q[0];
sx q[0];
rz(-3.0476373) q[0];
sx q[0];
rz(1.093338) q[0];
rz(0.023748485) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-2.645983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099305196) q[0];
sx q[0];
rz(-1.5047502) q[0];
sx q[0];
rz(1.5586067) q[0];
rz(-pi) q[1];
rz(0.42107196) q[2];
sx q[2];
rz(-2.6681528) q[2];
sx q[2];
rz(2.1630436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62140761) q[1];
sx q[1];
rz(-1.2222792) q[1];
sx q[1];
rz(-1.7640616) q[1];
rz(-2.7109259) q[3];
sx q[3];
rz(-1.1209295) q[3];
sx q[3];
rz(-1.3632185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89639837) q[2];
sx q[2];
rz(-0.79309016) q[2];
sx q[2];
rz(-2.8450361) q[2];
rz(0.5101997) q[3];
sx q[3];
rz(-1.6161796) q[3];
sx q[3];
rz(-0.27142522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940014) q[0];
sx q[0];
rz(-0.95518249) q[0];
sx q[0];
rz(1.1543132) q[0];
rz(1.1555903) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(-1.6824228) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302368) q[0];
sx q[0];
rz(-2.2855496) q[0];
sx q[0];
rz(-1.5404594) q[0];
x q[1];
rz(-2.295107) q[2];
sx q[2];
rz(-1.5716388) q[2];
sx q[2];
rz(-0.049132012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1016804) q[1];
sx q[1];
rz(-1.3915359) q[1];
sx q[1];
rz(1.0615361) q[1];
rz(-pi) q[2];
rz(1.0351712) q[3];
sx q[3];
rz(-1.0058306) q[3];
sx q[3];
rz(-1.9244827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3079754) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(-1.3758434) q[2];
rz(-0.086325072) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(-1.255792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1983222) q[0];
sx q[0];
rz(-2.4477796) q[0];
sx q[0];
rz(-0.86945239) q[0];
rz(0.72932875) q[1];
sx q[1];
rz(-2.2980233) q[1];
sx q[1];
rz(-2.9076911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9651523) q[0];
sx q[0];
rz(-1.5230706) q[0];
sx q[0];
rz(-0.41712572) q[0];
x q[1];
rz(-2.9314165) q[2];
sx q[2];
rz(-1.4769722) q[2];
sx q[2];
rz(1.891234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42185509) q[1];
sx q[1];
rz(-0.46793391) q[1];
sx q[1];
rz(0.33166064) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2083268) q[3];
sx q[3];
rz(-2.3911282) q[3];
sx q[3];
rz(-0.01736162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0539315) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.1762478) q[2];
rz(2.4704399) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(1.5161071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072356) q[0];
sx q[0];
rz(-0.86152995) q[0];
sx q[0];
rz(-1.0435411) q[0];
rz(0.097298233) q[1];
sx q[1];
rz(-2.0847335) q[1];
sx q[1];
rz(-1.1204488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1277043) q[0];
sx q[0];
rz(-2.5426425) q[0];
sx q[0];
rz(-1.6779792) q[0];
rz(-pi) q[1];
rz(0.78033041) q[2];
sx q[2];
rz(-2.188211) q[2];
sx q[2];
rz(1.2517901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35353684) q[1];
sx q[1];
rz(-0.18419838) q[1];
sx q[1];
rz(-2.8723529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9531526) q[3];
sx q[3];
rz(-1.095929) q[3];
sx q[3];
rz(-1.4120073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6843159) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(1.786001) q[2];
rz(1.5654303) q[3];
sx q[3];
rz(-2.0578945) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23715699) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(0.65658983) q[1];
sx q[1];
rz(-1.1142535) q[1];
sx q[1];
rz(2.5744892) q[1];
rz(2.6576853) q[2];
sx q[2];
rz(-2.172702) q[2];
sx q[2];
rz(-0.46401698) q[2];
rz(2.5593467) q[3];
sx q[3];
rz(-2.7979294) q[3];
sx q[3];
rz(2.9002849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
