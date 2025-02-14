OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(-2.1991576) q[1];
sx q[1];
rz(0.64185774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9489884) q[0];
sx q[0];
rz(-1.3697213) q[0];
sx q[0];
rz(-1.4549535) q[0];
rz(-pi) q[1];
rz(-1.5696758) q[2];
sx q[2];
rz(-1.5693451) q[2];
sx q[2];
rz(3.0641132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0986276) q[1];
sx q[1];
rz(-0.39448276) q[1];
sx q[1];
rz(-2.7261655) q[1];
rz(-pi) q[2];
rz(2.3083616) q[3];
sx q[3];
rz(-2.0249484) q[3];
sx q[3];
rz(-1.7851225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1610819) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(2.3316627) q[2];
rz(-0.03446456) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(-3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63475364) q[0];
sx q[0];
rz(-0.50951183) q[0];
sx q[0];
rz(2.4822045) q[0];
rz(-1.7022645) q[1];
sx q[1];
rz(-1.6292452) q[1];
sx q[1];
rz(-2.6606681) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.07671) q[0];
sx q[0];
rz(-1.1336898) q[0];
sx q[0];
rz(2.9124898) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3757295) q[2];
sx q[2];
rz(-0.89855902) q[2];
sx q[2];
rz(-1.6063251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2099328) q[1];
sx q[1];
rz(-1.5402964) q[1];
sx q[1];
rz(1.7790487) q[1];
x q[2];
rz(2.0388076) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(-0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7646358) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(2.5118206) q[2];
rz(-1.1791641) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(-2.0298957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(-2.1731398) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(-2.5586939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4752858) q[0];
sx q[0];
rz(-0.17284849) q[0];
sx q[0];
rz(-0.26474196) q[0];
rz(0.23704657) q[2];
sx q[2];
rz(-2.2893527) q[2];
sx q[2];
rz(1.3952554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0839094) q[1];
sx q[1];
rz(-2.1321802) q[1];
sx q[1];
rz(-0.76389216) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2887313) q[3];
sx q[3];
rz(-0.4133458) q[3];
sx q[3];
rz(1.0173544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53050238) q[2];
sx q[2];
rz(-3.0641596) q[2];
sx q[2];
rz(-0.21406847) q[2];
rz(0.30238447) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18836235) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(1.0444214) q[0];
rz(2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(2.9728319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8498189) q[0];
sx q[0];
rz(-2.5892604) q[0];
sx q[0];
rz(-2.0122347) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20511583) q[2];
sx q[2];
rz(-1.2502708) q[2];
sx q[2];
rz(1.6728624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9738439) q[1];
sx q[1];
rz(-0.97626057) q[1];
sx q[1];
rz(-1.3092625) q[1];
rz(-0.041958001) q[3];
sx q[3];
rz(-0.86554147) q[3];
sx q[3];
rz(0.77159568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3622482) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(-1.0603504) q[2];
rz(2.849071) q[3];
sx q[3];
rz(-2.4157603) q[3];
sx q[3];
rz(-0.084107548) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5523858) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(-2.2733083) q[0];
rz(-0.6262511) q[1];
sx q[1];
rz(-1.7528844) q[1];
sx q[1];
rz(1.7832322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5007965) q[0];
sx q[0];
rz(-1.8587041) q[0];
sx q[0];
rz(1.663371) q[0];
rz(-pi) q[1];
rz(-2.7600438) q[2];
sx q[2];
rz(-1.3693491) q[2];
sx q[2];
rz(-2.3214981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.061314) q[1];
sx q[1];
rz(-1.7931723) q[1];
sx q[1];
rz(2.910252) q[1];
rz(-2.0121947) q[3];
sx q[3];
rz(-0.98047653) q[3];
sx q[3];
rz(2.2143827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9354349) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(0.52331501) q[2];
rz(2.7715136) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.10941457) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(2.7874462) q[0];
rz(-0.94193637) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(-1.8249493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8982074) q[0];
sx q[0];
rz(-1.4665717) q[0];
sx q[0];
rz(2.6263155) q[0];
rz(-pi) q[1];
rz(-0.022158547) q[2];
sx q[2];
rz(-0.54931927) q[2];
sx q[2];
rz(-1.4469128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1209379) q[1];
sx q[1];
rz(-1.5096501) q[1];
sx q[1];
rz(0.17423363) q[1];
rz(-2.7735071) q[3];
sx q[3];
rz(-2.4637665) q[3];
sx q[3];
rz(-0.2673291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(2.6650688) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.7020096) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(1.4078377) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(2.5163311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1921034) q[0];
sx q[0];
rz(-1.2117627) q[0];
sx q[0];
rz(-0.60199317) q[0];
x q[1];
rz(-0.12282108) q[2];
sx q[2];
rz(-2.5993957) q[2];
sx q[2];
rz(1.2187097) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3846674) q[1];
sx q[1];
rz(-2.1672492) q[1];
sx q[1];
rz(1.3104964) q[1];
rz(-pi) q[2];
rz(-1.9373364) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(0.77223611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(1.7612877) q[2];
rz(-2.7050833) q[3];
sx q[3];
rz(-1.0218388) q[3];
sx q[3];
rz(2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9647144) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(0.51625133) q[0];
rz(2.3618354) q[1];
sx q[1];
rz(-0.98194352) q[1];
sx q[1];
rz(1.0468743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488432) q[0];
sx q[0];
rz(-1.7147281) q[0];
sx q[0];
rz(0.30152614) q[0];
x q[1];
rz(-0.5811695) q[2];
sx q[2];
rz(-2.0072674) q[2];
sx q[2];
rz(-1.798686) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2566764) q[1];
sx q[1];
rz(-0.11029989) q[1];
sx q[1];
rz(2.1259456) q[1];
rz(-pi) q[2];
rz(-0.67892142) q[3];
sx q[3];
rz(-0.58829868) q[3];
sx q[3];
rz(3.1117518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-2.011994) q[3];
sx q[3];
rz(-2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1421563) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(0.18705046) q[0];
rz(0.43141463) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(-1.6492708) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607912) q[0];
sx q[0];
rz(-1.6422762) q[0];
sx q[0];
rz(0.070822318) q[0];
x q[1];
rz(2.0772821) q[2];
sx q[2];
rz(-0.37831719) q[2];
sx q[2];
rz(1.0418721) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88384151) q[1];
sx q[1];
rz(-1.7353936) q[1];
sx q[1];
rz(1.1562892) q[1];
rz(-pi) q[2];
x q[2];
rz(3.047916) q[3];
sx q[3];
rz(-1.7266577) q[3];
sx q[3];
rz(0.81934281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6734068) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(-2.1569596) q[2];
rz(-2.6268688) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.2334067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.6615302) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(-1.1933391) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(-1.6903445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275996) q[0];
sx q[0];
rz(-1.7993449) q[0];
sx q[0];
rz(0.20968135) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7733712) q[2];
sx q[2];
rz(-1.612101) q[2];
sx q[2];
rz(-1.6172723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38498653) q[1];
sx q[1];
rz(-2.1861417) q[1];
sx q[1];
rz(0.032757515) q[1];
rz(-pi) q[2];
rz(-1.162446) q[3];
sx q[3];
rz(-1.2733885) q[3];
sx q[3];
rz(-1.4434467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5997368) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(0.80121458) q[2];
rz(2.2090705) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(0.41509375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(-0.32147944) q[1];
sx q[1];
rz(-0.97578661) q[1];
sx q[1];
rz(-1.6937561) q[1];
rz(2.9762815) q[2];
sx q[2];
rz(-1.6935692) q[2];
sx q[2];
rz(-1.6170681) q[2];
rz(0.23033167) q[3];
sx q[3];
rz(-1.2970222) q[3];
sx q[3];
rz(2.6877689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
