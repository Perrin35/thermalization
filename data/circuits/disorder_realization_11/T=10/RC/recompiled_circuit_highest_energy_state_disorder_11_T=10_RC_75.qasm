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
rz(-2.9304507) q[0];
sx q[0];
rz(-0.40786064) q[0];
rz(-1.7653699) q[1];
sx q[1];
rz(5.0566109) q[1];
sx q[1];
rz(15.772718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252547) q[0];
sx q[0];
rz(-3.0260575) q[0];
sx q[0];
rz(-1.6166685) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98528905) q[2];
sx q[2];
rz(-1.8409713) q[2];
sx q[2];
rz(0.5952685) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.237707) q[1];
sx q[1];
rz(-0.14169417) q[1];
sx q[1];
rz(-0.078448729) q[1];
x q[2];
rz(-0.14459855) q[3];
sx q[3];
rz(-2.5144092) q[3];
sx q[3];
rz(2.7759564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9593418) q[2];
sx q[2];
rz(-1.094341) q[2];
sx q[2];
rz(-2.0056637) q[2];
rz(-2.411339) q[3];
sx q[3];
rz(-1.7466931) q[3];
sx q[3];
rz(-1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092904329) q[0];
sx q[0];
rz(-0.95656675) q[0];
sx q[0];
rz(0.055305716) q[0];
rz(2.7703908) q[1];
sx q[1];
rz(-2.7822045) q[1];
sx q[1];
rz(-0.41195437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0028371) q[0];
sx q[0];
rz(-1.5056207) q[0];
sx q[0];
rz(0.52703339) q[0];
rz(-pi) q[1];
rz(1.2367593) q[2];
sx q[2];
rz(-1.8530117) q[2];
sx q[2];
rz(1.7441234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69393278) q[1];
sx q[1];
rz(-0.74697633) q[1];
sx q[1];
rz(-2.9332317) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92276908) q[3];
sx q[3];
rz(-1.2639746) q[3];
sx q[3];
rz(1.1359147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61341316) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(0.8826274) q[2];
rz(1.4334076) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(-0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9751137) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(2.574918) q[0];
rz(-0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(1.3272939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46337977) q[0];
sx q[0];
rz(-0.87317077) q[0];
sx q[0];
rz(0.11059983) q[0];
x q[1];
rz(-1.7812875) q[2];
sx q[2];
rz(-2.5062525) q[2];
sx q[2];
rz(2.6689305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9108563) q[1];
sx q[1];
rz(-1.5720782) q[1];
sx q[1];
rz(1.8788473) q[1];
rz(-2.684115) q[3];
sx q[3];
rz(-0.8972392) q[3];
sx q[3];
rz(-2.9411773) q[3];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4817151) q[0];
sx q[0];
rz(-0.91232038) q[0];
sx q[0];
rz(0.6657486) q[0];
rz(-2.2944229) q[1];
sx q[1];
rz(-0.24996346) q[1];
sx q[1];
rz(0.4096823) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9315791) q[0];
sx q[0];
rz(-1.0780452) q[0];
sx q[0];
rz(-1.3973622) q[0];
rz(0.3607765) q[2];
sx q[2];
rz(-0.84178001) q[2];
sx q[2];
rz(1.8052124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.332147) q[1];
sx q[1];
rz(-1.3464902) q[1];
sx q[1];
rz(-0.94384335) q[1];
rz(-pi) q[2];
x q[2];
rz(2.014854) q[3];
sx q[3];
rz(-1.5448625) q[3];
sx q[3];
rz(1.8166722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72982558) q[2];
sx q[2];
rz(-0.7080141) q[2];
sx q[2];
rz(1.7626308) q[2];
rz(2.0046558) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.1788504) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(0.41819292) q[0];
rz(-1.4232945) q[1];
sx q[1];
rz(-1.9530674) q[1];
sx q[1];
rz(1.6421912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2432775) q[0];
sx q[0];
rz(-2.0137798) q[0];
sx q[0];
rz(-1.8421296) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8078305) q[2];
sx q[2];
rz(-1.859772) q[2];
sx q[2];
rz(3.1372276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.027335701) q[1];
sx q[1];
rz(-1.789195) q[1];
sx q[1];
rz(0.059192358) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1058331) q[3];
sx q[3];
rz(-2.9966207) q[3];
sx q[3];
rz(1.5663869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2749918) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(-0.17821136) q[2];
rz(0.37402672) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(1.9988352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571534) q[0];
sx q[0];
rz(-1.8171808) q[0];
sx q[0];
rz(1.665218) q[0];
rz(2.5560675) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(-0.73385986) q[1];
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
x q[1];
rz(-1.8963054) q[2];
sx q[2];
rz(-1.2790888) q[2];
sx q[2];
rz(2.0355647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1230525) q[1];
sx q[1];
rz(-2.5830659) q[1];
sx q[1];
rz(-1.194242) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1185944) q[3];
sx q[3];
rz(-1.7590932) q[3];
sx q[3];
rz(-1.2533497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0840941) q[2];
sx q[2];
rz(-1.3749264) q[2];
sx q[2];
rz(1.6424087) q[2];
rz(-1.5205787) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(-2.9852941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(1.7479489) q[0];
rz(-2.5539894) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(0.29744068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2327691) q[0];
sx q[0];
rz(-0.89134514) q[0];
sx q[0];
rz(0.46436907) q[0];
x q[1];
rz(-1.6388551) q[2];
sx q[2];
rz(-0.27114332) q[2];
sx q[2];
rz(0.23672297) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1042034) q[1];
sx q[1];
rz(-0.89118499) q[1];
sx q[1];
rz(-0.012455539) q[1];
x q[2];
rz(2.2138811) q[3];
sx q[3];
rz(-1.3731908) q[3];
sx q[3];
rz(1.2090982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3191159) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(-2.759867) q[2];
rz(-0.60025674) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(-0.265358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(1.3375244) q[0];
rz(0.0094825347) q[1];
sx q[1];
rz(-2.1818706) q[1];
sx q[1];
rz(-1.3710075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66390304) q[0];
sx q[0];
rz(-0.36720095) q[0];
sx q[0];
rz(-0.3911256) q[0];
rz(-pi) q[1];
rz(2.1751733) q[2];
sx q[2];
rz(-0.76497173) q[2];
sx q[2];
rz(1.6993831) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1320033) q[1];
sx q[1];
rz(-1.4456962) q[1];
sx q[1];
rz(2.3522207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9389006) q[3];
sx q[3];
rz(-2.7198942) q[3];
sx q[3];
rz(-1.4055173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3235772) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(0.14249194) q[2];
rz(-2.2583708) q[3];
sx q[3];
rz(-1.3775237) q[3];
sx q[3];
rz(-1.6787329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.091081) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.2485414) q[0];
rz(-2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(-1.1936845) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4536205) q[0];
sx q[0];
rz(-1.2415452) q[0];
sx q[0];
rz(2.8696174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3825109) q[2];
sx q[2];
rz(-2.1164701) q[2];
sx q[2];
rz(-0.99198839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5387896) q[1];
sx q[1];
rz(-1.9059999) q[1];
sx q[1];
rz(0.523816) q[1];
rz(-pi) q[2];
rz(3.3832559e-05) q[3];
sx q[3];
rz(-1.3109929) q[3];
sx q[3];
rz(2.3855748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(-1.1129145) q[2];
rz(-3.0985966) q[3];
sx q[3];
rz(-1.6578511) q[3];
sx q[3];
rz(0.23785166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064903108) q[0];
sx q[0];
rz(-1.7273128) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(-1.9211357) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(-2.1641796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61924926) q[0];
sx q[0];
rz(-2.8641906) q[0];
sx q[0];
rz(2.7707556) q[0];
rz(1.8711817) q[2];
sx q[2];
rz(-1.4255878) q[2];
sx q[2];
rz(1.5086439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5557438) q[1];
sx q[1];
rz(-2.3341456) q[1];
sx q[1];
rz(2.2967352) q[1];
rz(0.39267003) q[3];
sx q[3];
rz(-1.8373946) q[3];
sx q[3];
rz(-0.50470282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0105373) q[2];
sx q[2];
rz(-2.8752893) q[2];
sx q[2];
rz(-1.6923426) q[2];
rz(0.87861711) q[3];
sx q[3];
rz(-2.0802458) q[3];
sx q[3];
rz(1.3531468) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286248) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.3799008) q[1];
sx q[1];
rz(-1.6086144) q[1];
sx q[1];
rz(-0.10133941) q[1];
rz(2.9607282) q[2];
sx q[2];
rz(-1.3915863) q[2];
sx q[2];
rz(-2.1044921) q[2];
rz(-2.8158549) q[3];
sx q[3];
rz(-2.4690463) q[3];
sx q[3];
rz(1.9161968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
