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
rz(2.9003484) q[0];
sx q[0];
rz(3.3527346) q[0];
sx q[0];
rz(9.0169173) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(-0.064755138) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714343) q[0];
sx q[0];
rz(-1.4553833) q[0];
sx q[0];
rz(0.0053216318) q[0];
rz(-2.8207755) q[2];
sx q[2];
rz(-1.0091558) q[2];
sx q[2];
rz(2.3412242) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3169466) q[1];
sx q[1];
rz(-1.4295409) q[1];
sx q[1];
rz(-1.5596175) q[1];
x q[2];
rz(2.519386) q[3];
sx q[3];
rz(-1.4861306) q[3];
sx q[3];
rz(2.0537927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18225081) q[2];
sx q[2];
rz(-1.094341) q[2];
sx q[2];
rz(-1.135929) q[2];
rz(2.411339) q[3];
sx q[3];
rz(-1.7466931) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092904329) q[0];
sx q[0];
rz(-2.1850259) q[0];
sx q[0];
rz(-3.0862869) q[0];
rz(-0.37120184) q[1];
sx q[1];
rz(-2.7822045) q[1];
sx q[1];
rz(-0.41195437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7474351) q[0];
sx q[0];
rz(-2.0965946) q[0];
sx q[0];
rz(-1.4954241) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2949176) q[2];
sx q[2];
rz(-2.7077423) q[2];
sx q[2];
rz(-2.2920319) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69393278) q[1];
sx q[1];
rz(-2.3946163) q[1];
sx q[1];
rz(-2.9332317) q[1];
x q[2];
rz(-2.2188236) q[3];
sx q[3];
rz(-1.8776181) q[3];
sx q[3];
rz(1.1359147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61341316) q[2];
sx q[2];
rz(-1.6772905) q[2];
sx q[2];
rz(-0.8826274) q[2];
rz(-1.708185) q[3];
sx q[3];
rz(-1.9288758) q[3];
sx q[3];
rz(-2.9679756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-0.16647896) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(0.56667462) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.5387225) q[1];
sx q[1];
rz(1.3272939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2921995) q[0];
sx q[0];
rz(-0.70488323) q[0];
sx q[0];
rz(-1.439875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3603052) q[2];
sx q[2];
rz(-2.5062525) q[2];
sx q[2];
rz(2.6689305) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23073634) q[1];
sx q[1];
rz(-1.5695144) q[1];
sx q[1];
rz(1.2627453) q[1];
rz(-0.8437952) q[3];
sx q[3];
rz(-1.9232755) q[3];
sx q[3];
rz(-1.6683287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3920445) q[2];
sx q[2];
rz(-2.0508524) q[2];
sx q[2];
rz(-2.2331179) q[2];
rz(0.79489094) q[3];
sx q[3];
rz(-1.1029693) q[3];
sx q[3];
rz(-3.095043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817151) q[0];
sx q[0];
rz(-2.2292723) q[0];
sx q[0];
rz(-2.4758441) q[0];
rz(2.2944229) q[1];
sx q[1];
rz(-2.8916292) q[1];
sx q[1];
rz(-2.7319103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6981229) q[0];
sx q[0];
rz(-1.4181678) q[0];
sx q[0];
rz(2.6425155) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7808162) q[2];
sx q[2];
rz(-2.2998126) q[2];
sx q[2];
rz(1.3363802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.078881513) q[1];
sx q[1];
rz(-0.96187519) q[1];
sx q[1];
rz(-0.2746065) q[1];
x q[2];
rz(3.112875) q[3];
sx q[3];
rz(-2.014694) q[3];
sx q[3];
rz(2.9080528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4117671) q[2];
sx q[2];
rz(-0.7080141) q[2];
sx q[2];
rz(-1.3789619) q[2];
rz(1.1369368) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(-0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9627422) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(-2.7233997) q[0];
rz(-1.4232945) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(-1.6421912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8189637) q[0];
sx q[0];
rz(-0.51473619) q[0];
sx q[0];
rz(-2.6273652) q[0];
rz(1.8756494) q[2];
sx q[2];
rz(-1.2513759) q[2];
sx q[2];
rz(-1.4679421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5852927) q[1];
sx q[1];
rz(-1.628581) q[1];
sx q[1];
rz(-1.3520266) q[1];
rz(1.70056) q[3];
sx q[3];
rz(-1.6356182) q[3];
sx q[3];
rz(-2.685252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8666009) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(2.9633813) q[2];
rz(0.37402672) q[3];
sx q[3];
rz(-0.97291294) q[3];
sx q[3];
rz(-1.9988352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571534) q[0];
sx q[0];
rz(-1.3244119) q[0];
sx q[0];
rz(-1.665218) q[0];
rz(-0.58552512) q[1];
sx q[1];
rz(-2.245677) q[1];
sx q[1];
rz(-2.4077328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.748679) q[0];
sx q[0];
rz(-0.80930149) q[0];
sx q[0];
rz(0.63776638) q[0];
rz(2.8346905) q[2];
sx q[2];
rz(-1.2595121) q[2];
sx q[2];
rz(-0.5615304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3702995) q[1];
sx q[1];
rz(-1.3746737) q[1];
sx q[1];
rz(2.0972154) q[1];
x q[2];
rz(1.9820342) q[3];
sx q[3];
rz(-0.48732584) q[3];
sx q[3];
rz(-3.0912378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.057498589) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(1.499184) q[2];
rz(-1.621014) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(-0.15629855) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(1.7479489) q[0];
rz(-0.58760324) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(2.844152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2327691) q[0];
sx q[0];
rz(-2.2502475) q[0];
sx q[0];
rz(0.46436907) q[0];
rz(-pi) q[1];
rz(-1.6388551) q[2];
sx q[2];
rz(-2.8704493) q[2];
sx q[2];
rz(-0.23672297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0373893) q[1];
sx q[1];
rz(-2.2504077) q[1];
sx q[1];
rz(3.1291371) q[1];
x q[2];
rz(-1.8930412) q[3];
sx q[3];
rz(-0.6686223) q[3];
sx q[3];
rz(-3.0361255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3191159) q[2];
sx q[2];
rz(-1.7812704) q[2];
sx q[2];
rz(2.759867) q[2];
rz(-2.5413359) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(-2.8762347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47939077) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(1.3375244) q[0];
rz(-0.0094825347) q[1];
sx q[1];
rz(-2.1818706) q[1];
sx q[1];
rz(-1.7705852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6021045) q[0];
sx q[0];
rz(-1.7080902) q[0];
sx q[0];
rz(2.7999393) q[0];
x q[1];
rz(-0.49937906) q[2];
sx q[2];
rz(-0.96448318) q[2];
sx q[2];
rz(0.93580907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43804719) q[1];
sx q[1];
rz(-2.3444972) q[1];
sx q[1];
rz(-2.966267) q[1];
x q[2];
rz(-1.6608606) q[3];
sx q[3];
rz(-1.9833296) q[3];
sx q[3];
rz(1.6270669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3235772) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(2.9991007) q[2];
rz(2.2583708) q[3];
sx q[3];
rz(-1.3775237) q[3];
sx q[3];
rz(-1.4628598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0505117) q[0];
sx q[0];
rz(-3.0497146) q[0];
sx q[0];
rz(1.8930513) q[0];
rz(-2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(-1.1936845) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4536205) q[0];
sx q[0];
rz(-1.2415452) q[0];
sx q[0];
rz(-0.27197522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72288469) q[2];
sx q[2];
rz(-0.9019081) q[2];
sx q[2];
rz(-1.0793874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.602803) q[1];
sx q[1];
rz(-1.9059999) q[1];
sx q[1];
rz(2.6177767) q[1];
rz(-pi) q[2];
rz(1.570669) q[3];
sx q[3];
rz(-2.8817892) q[3];
sx q[3];
rz(0.75614959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(2.0286782) q[2];
rz(3.0985966) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-2.903741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766895) q[0];
sx q[0];
rz(-1.4142798) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(-1.2204569) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(-0.97741309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138151) q[0];
sx q[0];
rz(-1.8288946) q[0];
sx q[0];
rz(-1.6736223) q[0];
rz(-pi) q[1];
rz(1.270411) q[2];
sx q[2];
rz(-1.4255878) q[2];
sx q[2];
rz(1.6329488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7070476) q[1];
sx q[1];
rz(-1.0705528) q[1];
sx q[1];
rz(0.90737271) q[1];
x q[2];
rz(2.5217358) q[3];
sx q[3];
rz(-2.6708948) q[3];
sx q[3];
rz(2.6420267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1310554) q[2];
sx q[2];
rz(-2.8752893) q[2];
sx q[2];
rz(1.6923426) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-1.0613469) q[3];
sx q[3];
rz(1.3531468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286248) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.7616918) q[1];
sx q[1];
rz(-1.5329783) q[1];
sx q[1];
rz(3.0402532) q[1];
rz(-2.352667) q[2];
sx q[2];
rz(-2.8876705) q[2];
sx q[2];
rz(-1.3063274) q[2];
rz(-2.4951999) q[3];
sx q[3];
rz(-1.7714995) q[3];
sx q[3];
rz(-2.5378791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
