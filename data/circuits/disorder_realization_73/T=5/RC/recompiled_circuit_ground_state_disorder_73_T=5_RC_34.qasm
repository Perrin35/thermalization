OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1239399) q[0];
sx q[0];
rz(-0.91683638) q[0];
sx q[0];
rz(-0.43489781) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(3.6606001) q[1];
sx q[1];
rz(11.984348) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911444) q[0];
sx q[0];
rz(-1.309937) q[0];
sx q[0];
rz(2.9820739) q[0];
rz(-0.93272722) q[2];
sx q[2];
rz(-1.9550394) q[2];
sx q[2];
rz(0.24894938) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4705321) q[1];
sx q[1];
rz(-1.8166421) q[1];
sx q[1];
rz(-1.125294) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1277027) q[3];
sx q[3];
rz(-2.8134973) q[3];
sx q[3];
rz(1.2893639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55328289) q[2];
sx q[2];
rz(-1.8841691) q[2];
sx q[2];
rz(-0.29169875) q[2];
rz(2.6907673) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(0.60744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4400103) q[0];
sx q[0];
rz(-0.99047438) q[0];
sx q[0];
rz(-2.0462346) q[0];
rz(1.9502684) q[1];
sx q[1];
rz(-0.92000633) q[1];
sx q[1];
rz(-2.2390168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94287017) q[0];
sx q[0];
rz(-0.96317055) q[0];
sx q[0];
rz(-0.81815079) q[0];
rz(-0.90229374) q[2];
sx q[2];
rz(-1.3375139) q[2];
sx q[2];
rz(1.4177314) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7831259) q[1];
sx q[1];
rz(-1.2791113) q[1];
sx q[1];
rz(-1.6963708) q[1];
x q[2];
rz(0.9489305) q[3];
sx q[3];
rz(-1.9919259) q[3];
sx q[3];
rz(-1.7647183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9537182) q[2];
sx q[2];
rz(-2.1672921) q[2];
sx q[2];
rz(3.0226959) q[2];
rz(2.4925354) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.4547179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4897937) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(-0.5157665) q[0];
rz(2.0491397) q[1];
sx q[1];
rz(-0.40320435) q[1];
sx q[1];
rz(1.8663503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1920722) q[0];
sx q[0];
rz(-1.6123591) q[0];
sx q[0];
rz(-1.5526718) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3408717) q[2];
sx q[2];
rz(-1.2650239) q[2];
sx q[2];
rz(1.1912322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67120463) q[1];
sx q[1];
rz(-0.71031308) q[1];
sx q[1];
rz(1.2673668) q[1];
rz(0.79408349) q[3];
sx q[3];
rz(-1.6424254) q[3];
sx q[3];
rz(-0.41251999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.065993) q[2];
sx q[2];
rz(-1.8241355) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(-0.48948151) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(2.421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2751665) q[0];
sx q[0];
rz(-0.32298276) q[0];
sx q[0];
rz(0.48165709) q[0];
rz(0.9306759) q[1];
sx q[1];
rz(-1.296867) q[1];
sx q[1];
rz(3.1226588) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16068072) q[0];
sx q[0];
rz(-1.5917516) q[0];
sx q[0];
rz(1.6208072) q[0];
rz(1.5140947) q[2];
sx q[2];
rz(-2.4783274) q[2];
sx q[2];
rz(-0.5723638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1065797) q[1];
sx q[1];
rz(-2.0845319) q[1];
sx q[1];
rz(0.68961838) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75441672) q[3];
sx q[3];
rz(-2.2212611) q[3];
sx q[3];
rz(0.90620302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21861741) q[2];
sx q[2];
rz(-2.7851892) q[2];
sx q[2];
rz(2.3962928) q[2];
rz(-2.7636012) q[3];
sx q[3];
rz(-1.1443006) q[3];
sx q[3];
rz(2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22076631) q[0];
sx q[0];
rz(-0.31196088) q[0];
sx q[0];
rz(1.1908603) q[0];
rz(-2.6530755) q[1];
sx q[1];
rz(-2.2256475) q[1];
sx q[1];
rz(0.8078422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69881781) q[0];
sx q[0];
rz(-1.6614058) q[0];
sx q[0];
rz(-1.242799) q[0];
x q[1];
rz(-0.22831445) q[2];
sx q[2];
rz(-2.0234152) q[2];
sx q[2];
rz(1.7646947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74238363) q[1];
sx q[1];
rz(-2.0258372) q[1];
sx q[1];
rz(-1.6687523) q[1];
x q[2];
rz(0.2000974) q[3];
sx q[3];
rz(-2.1842624) q[3];
sx q[3];
rz(0.6544978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0025582) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(-2.9774418) q[2];
rz(-0.74603355) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(0.073808864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98110759) q[0];
sx q[0];
rz(-1.8758513) q[0];
sx q[0];
rz(2.4142081) q[0];
rz(0.47850594) q[1];
sx q[1];
rz(-1.452927) q[1];
sx q[1];
rz(-0.7368288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751303) q[0];
sx q[0];
rz(-1.7084165) q[0];
sx q[0];
rz(3.0855975) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1428035) q[2];
sx q[2];
rz(-0.30747947) q[2];
sx q[2];
rz(0.82915074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7718862) q[1];
sx q[1];
rz(-0.73644887) q[1];
sx q[1];
rz(2.1320029) q[1];
rz(-pi) q[2];
x q[2];
rz(0.013929587) q[3];
sx q[3];
rz(-1.9584112) q[3];
sx q[3];
rz(-2.882021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9937146) q[2];
sx q[2];
rz(-1.0241877) q[2];
sx q[2];
rz(-2.4564339) q[2];
rz(2.1327175) q[3];
sx q[3];
rz(-0.69005552) q[3];
sx q[3];
rz(0.25079295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93609) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(1.093338) q[0];
rz(-0.023748485) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-0.49560961) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4706866) q[0];
sx q[0];
rz(-1.5829594) q[0];
sx q[0];
rz(-3.0755416) q[0];
x q[1];
rz(2.7042386) q[2];
sx q[2];
rz(-1.7582571) q[2];
sx q[2];
rz(2.1700493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1417947) q[1];
sx q[1];
rz(-0.3965946) q[1];
sx q[1];
rz(-0.48626089) q[1];
x q[2];
rz(-1.0823332) q[3];
sx q[3];
rz(-1.1853855) q[3];
sx q[3];
rz(3.1311991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(0.29655656) q[2];
rz(2.631393) q[3];
sx q[3];
rz(-1.6161796) q[3];
sx q[3];
rz(-2.8701674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940014) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(1.1543132) q[0];
rz(-1.1555903) q[1];
sx q[1];
rz(-0.65735936) q[1];
sx q[1];
rz(1.4591699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66508913) q[0];
sx q[0];
rz(-0.71528331) q[0];
sx q[0];
rz(0.034937783) q[0];
rz(-pi) q[1];
rz(-3.1404678) q[2];
sx q[2];
rz(-2.2951067) q[2];
sx q[2];
rz(1.520919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0399122) q[1];
sx q[1];
rz(-1.3915359) q[1];
sx q[1];
rz(2.0800566) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67787637) q[3];
sx q[3];
rz(-0.75784412) q[3];
sx q[3];
rz(2.0539396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83361721) q[2];
sx q[2];
rz(-2.5994382) q[2];
sx q[2];
rz(-1.3758434) q[2];
rz(3.0552676) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(1.8858006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9432705) q[0];
sx q[0];
rz(-0.69381303) q[0];
sx q[0];
rz(2.2721403) q[0];
rz(-0.72932875) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(0.23390153) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37321842) q[0];
sx q[0];
rz(-1.1541751) q[0];
sx q[0];
rz(1.6229902) q[0];
rz(-pi) q[1];
rz(-0.42371427) q[2];
sx q[2];
rz(-0.22988453) q[2];
sx q[2];
rz(2.4073441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7197376) q[1];
sx q[1];
rz(-0.46793391) q[1];
sx q[1];
rz(-2.809932) q[1];
rz(0.31932217) q[3];
sx q[3];
rz(-0.87933137) q[3];
sx q[3];
rz(-2.6806074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.087661155) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(-1.9653448) q[2];
rz(2.4704399) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(1.5161071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072356) q[0];
sx q[0];
rz(-2.2800627) q[0];
sx q[0];
rz(1.0435411) q[0];
rz(3.0442944) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(2.0211438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0138884) q[0];
sx q[0];
rz(-2.5426425) q[0];
sx q[0];
rz(1.6779792) q[0];
x q[1];
rz(-0.78033041) q[2];
sx q[2];
rz(-0.95338168) q[2];
sx q[2];
rz(1.2517901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6594362) q[1];
sx q[1];
rz(-1.6195356) q[1];
sx q[1];
rz(-2.963889) q[1];
rz(-pi) q[2];
rz(-2.9531526) q[3];
sx q[3];
rz(-2.0456637) q[3];
sx q[3];
rz(1.7295854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6843159) q[2];
sx q[2];
rz(-0.61351675) q[2];
sx q[2];
rz(-1.3555917) q[2];
rz(1.5654303) q[3];
sx q[3];
rz(-2.0578945) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044357) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(2.4850028) q[1];
sx q[1];
rz(-2.0273392) q[1];
sx q[1];
rz(-0.56710342) q[1];
rz(0.48390735) q[2];
sx q[2];
rz(-0.96889062) q[2];
sx q[2];
rz(2.6775757) q[2];
rz(-2.5593467) q[3];
sx q[3];
rz(-0.34366321) q[3];
sx q[3];
rz(-0.2413078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
