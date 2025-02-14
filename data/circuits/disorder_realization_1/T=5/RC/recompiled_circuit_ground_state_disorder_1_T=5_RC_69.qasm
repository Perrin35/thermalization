OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.863707) q[0];
sx q[0];
rz(-0.42930332) q[0];
sx q[0];
rz(0.31097809) q[0];
rz(1.7929945) q[1];
sx q[1];
rz(-1.0962948) q[1];
sx q[1];
rz(0.57408339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047217) q[0];
sx q[0];
rz(-1.939491) q[0];
sx q[0];
rz(1.1218907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1537504) q[2];
sx q[2];
rz(-1.541809) q[2];
sx q[2];
rz(1.0430973) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8771389) q[1];
sx q[1];
rz(-1.2065556) q[1];
sx q[1];
rz(-3.0458782) q[1];
x q[2];
rz(-1.9472136) q[3];
sx q[3];
rz(-1.6302846) q[3];
sx q[3];
rz(0.74947442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2744039) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(1.8223507) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.8354445) q[3];
sx q[3];
rz(-2.3225972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(2.4617885) q[0];
rz(-0.80348429) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(0.63978535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6006192) q[0];
sx q[0];
rz(-0.94087155) q[0];
sx q[0];
rz(2.7472904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4117208) q[2];
sx q[2];
rz(-1.4966855) q[2];
sx q[2];
rz(-0.48219901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3099932) q[1];
sx q[1];
rz(-0.91504708) q[1];
sx q[1];
rz(2.3393231) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62044797) q[3];
sx q[3];
rz(-2.7683966) q[3];
sx q[3];
rz(-0.20388145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60376072) q[2];
sx q[2];
rz(-0.49850264) q[2];
sx q[2];
rz(-2.4375088) q[2];
rz(3.0294561) q[3];
sx q[3];
rz(-1.4487368) q[3];
sx q[3];
rz(2.1191547) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10244399) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(2.802134) q[0];
rz(-2.0231694) q[1];
sx q[1];
rz(-0.92670852) q[1];
sx q[1];
rz(3.1059473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0775722) q[0];
sx q[0];
rz(-1.7776995) q[0];
sx q[0];
rz(-1.6585361) q[0];
rz(-pi) q[1];
rz(1.8730803) q[2];
sx q[2];
rz(-1.8955823) q[2];
sx q[2];
rz(-3.0168282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.798637) q[1];
sx q[1];
rz(-1.7854714) q[1];
sx q[1];
rz(0.66185419) q[1];
x q[2];
rz(0.6142637) q[3];
sx q[3];
rz(-1.302056) q[3];
sx q[3];
rz(2.4128071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79225916) q[2];
sx q[2];
rz(-0.68816853) q[2];
sx q[2];
rz(0.86000693) q[2];
rz(2.3490014) q[3];
sx q[3];
rz(-2.9383797) q[3];
sx q[3];
rz(2.121076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0122796) q[0];
sx q[0];
rz(-1.3849994) q[0];
sx q[0];
rz(-0.63283515) q[0];
rz(-1.9889779) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(1.4702183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7790079) q[0];
sx q[0];
rz(-2.2200395) q[0];
sx q[0];
rz(-3.1253689) q[0];
x q[1];
rz(-0.15369065) q[2];
sx q[2];
rz(-2.3665303) q[2];
sx q[2];
rz(-2.3714921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1296334) q[1];
sx q[1];
rz(-2.2599038) q[1];
sx q[1];
rz(-0.33504519) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4864462) q[3];
sx q[3];
rz(-1.0628502) q[3];
sx q[3];
rz(-0.20482132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70964885) q[2];
sx q[2];
rz(-0.96082965) q[2];
sx q[2];
rz(2.7107837) q[2];
rz(2.4591947) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(-0.94990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34669852) q[0];
sx q[0];
rz(-1.2656724) q[0];
sx q[0];
rz(1.9516113) q[0];
rz(2.8311912) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(0.50500542) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0895665) q[0];
sx q[0];
rz(-1.0086233) q[0];
sx q[0];
rz(-2.6181391) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0104899) q[2];
sx q[2];
rz(-1.982568) q[2];
sx q[2];
rz(-2.6302261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47520275) q[1];
sx q[1];
rz(-1.0745) q[1];
sx q[1];
rz(1.6700891) q[1];
rz(-0.72568958) q[3];
sx q[3];
rz(-1.5066875) q[3];
sx q[3];
rz(-0.025246092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7348822) q[2];
sx q[2];
rz(-0.70047417) q[2];
sx q[2];
rz(-2.238838) q[2];
rz(-2.086575) q[3];
sx q[3];
rz(-0.86692923) q[3];
sx q[3];
rz(-2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91586739) q[0];
sx q[0];
rz(-2.4287455) q[0];
sx q[0];
rz(2.0676887) q[0];
rz(-1.7621) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(0.097537907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65501761) q[0];
sx q[0];
rz(-1.7518105) q[0];
sx q[0];
rz(-0.54365309) q[0];
x q[1];
rz(-2.9729207) q[2];
sx q[2];
rz(-1.198581) q[2];
sx q[2];
rz(1.3160567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8360459) q[1];
sx q[1];
rz(-1.6191481) q[1];
sx q[1];
rz(-1.509686) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.318368) q[3];
sx q[3];
rz(-0.71079094) q[3];
sx q[3];
rz(1.2430134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11233106) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(-2.7304999) q[2];
rz(-1.7279651) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(-1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.3635062) q[0];
sx q[0];
rz(-3.111105) q[0];
sx q[0];
rz(1.4084858) q[0];
rz(-2.1259437) q[1];
sx q[1];
rz(-1.8135095) q[1];
sx q[1];
rz(1.6633063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7218839) q[0];
sx q[0];
rz(-0.96707487) q[0];
sx q[0];
rz(-1.9323856) q[0];
rz(-pi) q[1];
rz(2.4934216) q[2];
sx q[2];
rz(-1.0888466) q[2];
sx q[2];
rz(2.3881557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1455545) q[1];
sx q[1];
rz(-2.013447) q[1];
sx q[1];
rz(0.084060479) q[1];
x q[2];
rz(-0.13161259) q[3];
sx q[3];
rz(-2.0672879) q[3];
sx q[3];
rz(-2.541996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0279072) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(0.17024635) q[2];
rz(-1.8611192) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(1.1580275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-1.4729507) q[0];
sx q[0];
rz(-0.96422115) q[0];
sx q[0];
rz(-2.6131795) q[0];
rz(0.93327418) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(-2.5859213) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816102) q[0];
sx q[0];
rz(-2.1851808) q[0];
sx q[0];
rz(-0.23991628) q[0];
rz(-pi) q[1];
rz(1.6711177) q[2];
sx q[2];
rz(-2.360376) q[2];
sx q[2];
rz(0.90512102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5720776) q[1];
sx q[1];
rz(-1.5387427) q[1];
sx q[1];
rz(2.4299014) q[1];
x q[2];
rz(-1.2597841) q[3];
sx q[3];
rz(-0.90299435) q[3];
sx q[3];
rz(0.0086812191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8133424) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(-2.8343406) q[2];
rz(0.79636374) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55427134) q[0];
sx q[0];
rz(-2.0926496) q[0];
sx q[0];
rz(1.8121423) q[0];
rz(-2.9216271) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(-1.3185917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466965) q[0];
sx q[0];
rz(-1.7864972) q[0];
sx q[0];
rz(1.8442276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21130113) q[2];
sx q[2];
rz(-0.50926477) q[2];
sx q[2];
rz(-2.3788521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4107549) q[1];
sx q[1];
rz(-1.9111834) q[1];
sx q[1];
rz(-2.0654668) q[1];
rz(-pi) q[2];
rz(0.52826442) q[3];
sx q[3];
rz(-2.3481728) q[3];
sx q[3];
rz(-1.8509282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69507504) q[2];
sx q[2];
rz(-1.4277642) q[2];
sx q[2];
rz(-1.8565149) q[2];
rz(-0.88589969) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.6133962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0415683) q[0];
sx q[0];
rz(-0.73606857) q[0];
sx q[0];
rz(-2.8247483) q[0];
rz(-2.7007804) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(-2.7712834) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99890868) q[0];
sx q[0];
rz(-0.38247358) q[0];
sx q[0];
rz(-1.6539025) q[0];
rz(-pi) q[1];
rz(-2.2971588) q[2];
sx q[2];
rz(-2.0877247) q[2];
sx q[2];
rz(-0.58973613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1052484) q[1];
sx q[1];
rz(-0.40648983) q[1];
sx q[1];
rz(0.84291934) q[1];
x q[2];
rz(-1.466981) q[3];
sx q[3];
rz(-1.6214288) q[3];
sx q[3];
rz(1.6679483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4050498) q[2];
sx q[2];
rz(-0.49095792) q[2];
sx q[2];
rz(1.5092108) q[2];
rz(2.3368808) q[3];
sx q[3];
rz(-1.6377621) q[3];
sx q[3];
rz(-0.10778431) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1152773) q[0];
sx q[0];
rz(-1.7290709) q[0];
sx q[0];
rz(1.2363634) q[0];
rz(0.23030494) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(-0.7028107) q[2];
sx q[2];
rz(-2.3431449) q[2];
sx q[2];
rz(0.14771067) q[2];
rz(-1.1272011) q[3];
sx q[3];
rz(-2.4738635) q[3];
sx q[3];
rz(0.86480468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
