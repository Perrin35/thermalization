OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2603962) q[0];
sx q[0];
rz(1.908778) q[0];
sx q[0];
rz(9.0969152) q[0];
rz(2.7479563) q[1];
sx q[1];
rz(-1.6180232) q[1];
sx q[1];
rz(1.5394428) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6820591) q[0];
sx q[0];
rz(-1.5757808) q[0];
sx q[0];
rz(2.940516) q[0];
x q[1];
rz(1.1423768) q[2];
sx q[2];
rz(-0.5912458) q[2];
sx q[2];
rz(-0.57280409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0346759) q[1];
sx q[1];
rz(-2.6991416) q[1];
sx q[1];
rz(-0.4974858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39253916) q[3];
sx q[3];
rz(-1.919636) q[3];
sx q[3];
rz(-0.67482812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9615122) q[2];
sx q[2];
rz(-2.1891258) q[2];
sx q[2];
rz(-2.3251779) q[2];
rz(-0.69711971) q[3];
sx q[3];
rz(-2.7643047) q[3];
sx q[3];
rz(-2.2009946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546251) q[0];
sx q[0];
rz(-2.0836199) q[0];
sx q[0];
rz(0.70607287) q[0];
rz(0.1708897) q[1];
sx q[1];
rz(-0.51097521) q[1];
sx q[1];
rz(-2.1461646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0577431) q[0];
sx q[0];
rz(-0.47561663) q[0];
sx q[0];
rz(2.3887079) q[0];
rz(-pi) q[1];
rz(-2.1336251) q[2];
sx q[2];
rz(-0.78608222) q[2];
sx q[2];
rz(-0.85725923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.625841) q[1];
sx q[1];
rz(-1.4179572) q[1];
sx q[1];
rz(-1.3550853) q[1];
x q[2];
rz(0.63748116) q[3];
sx q[3];
rz(-1.5439171) q[3];
sx q[3];
rz(2.0826552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.55676111) q[2];
sx q[2];
rz(-0.91412592) q[2];
sx q[2];
rz(-2.8625028) q[2];
rz(-2.4537405) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(-2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2897709) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(-2.7626792) q[0];
rz(-1.9508349) q[1];
sx q[1];
rz(-2.7749116) q[1];
sx q[1];
rz(0.77622882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.533894) q[0];
sx q[0];
rz(-1.726681) q[0];
sx q[0];
rz(1.1381696) q[0];
rz(-pi) q[1];
rz(1.4188781) q[2];
sx q[2];
rz(-1.3891274) q[2];
sx q[2];
rz(2.7318397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8175428) q[1];
sx q[1];
rz(-2.2607798) q[1];
sx q[1];
rz(1.2381366) q[1];
rz(-pi) q[2];
rz(-1.0465195) q[3];
sx q[3];
rz(-1.5308497) q[3];
sx q[3];
rz(-2.0986882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1227405) q[2];
sx q[2];
rz(-1.9032225) q[2];
sx q[2];
rz(-2.5416601) q[2];
rz(-1.2395202) q[3];
sx q[3];
rz(-1.4116838) q[3];
sx q[3];
rz(0.6111353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83830738) q[0];
sx q[0];
rz(-2.4628283) q[0];
sx q[0];
rz(-1.7093866) q[0];
rz(2.2749061) q[1];
sx q[1];
rz(-1.3431834) q[1];
sx q[1];
rz(1.0923045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9878532) q[0];
sx q[0];
rz(-0.9747552) q[0];
sx q[0];
rz(-2.2558262) q[0];
x q[1];
rz(2.1020911) q[2];
sx q[2];
rz(-2.4091125) q[2];
sx q[2];
rz(-1.3074444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0835613) q[1];
sx q[1];
rz(-1.5688854) q[1];
sx q[1];
rz(0.46431904) q[1];
rz(-2.3471429) q[3];
sx q[3];
rz(-2.2243119) q[3];
sx q[3];
rz(0.78032035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0563353) q[2];
sx q[2];
rz(-1.9838355) q[2];
sx q[2];
rz(1.2752656) q[2];
rz(1.3750252) q[3];
sx q[3];
rz(-1.9019889) q[3];
sx q[3];
rz(2.1084771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035868693) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(1.3637654) q[0];
rz(-2.5719602) q[1];
sx q[1];
rz(-2.6424347) q[1];
sx q[1];
rz(2.6380576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0346876) q[0];
sx q[0];
rz(-1.4464753) q[0];
sx q[0];
rz(2.9479545) q[0];
rz(-pi) q[1];
rz(2.4120008) q[2];
sx q[2];
rz(-0.9895095) q[2];
sx q[2];
rz(0.94129291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8263957) q[1];
sx q[1];
rz(-0.082210899) q[1];
sx q[1];
rz(-1.7226687) q[1];
x q[2];
rz(0.77025099) q[3];
sx q[3];
rz(-1.3512843) q[3];
sx q[3];
rz(-1.3971922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86294952) q[2];
sx q[2];
rz(-2.7549665) q[2];
sx q[2];
rz(-1.6031727) q[2];
rz(3.0053511) q[3];
sx q[3];
rz(-1.7377661) q[3];
sx q[3];
rz(-1.9041825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.49029008) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(-0.34143099) q[0];
rz(0.731172) q[1];
sx q[1];
rz(-0.57594222) q[1];
sx q[1];
rz(2.0700571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51136953) q[0];
sx q[0];
rz(-1.9804738) q[0];
sx q[0];
rz(1.131534) q[0];
rz(-pi) q[1];
rz(-3.0546636) q[2];
sx q[2];
rz(-2.5247716) q[2];
sx q[2];
rz(1.4035743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.092128828) q[1];
sx q[1];
rz(-0.84869203) q[1];
sx q[1];
rz(-2.5404088) q[1];
x q[2];
rz(-1.5626181) q[3];
sx q[3];
rz(-2.6371752) q[3];
sx q[3];
rz(0.8917419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3557446) q[2];
sx q[2];
rz(-1.637849) q[2];
sx q[2];
rz(-0.49622932) q[2];
rz(-1.5931386) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(1.0065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9736495) q[0];
sx q[0];
rz(-1.668546) q[0];
sx q[0];
rz(3.1267401) q[0];
rz(-1.4133833) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(-0.76494876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0293065) q[0];
sx q[0];
rz(-0.95369442) q[0];
sx q[0];
rz(-0.86227148) q[0];
rz(-pi) q[1];
rz(1.2194446) q[2];
sx q[2];
rz(-0.91809154) q[2];
sx q[2];
rz(-1.2437133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5523259) q[1];
sx q[1];
rz(-1.5752042) q[1];
sx q[1];
rz(0.39673504) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9243605) q[3];
sx q[3];
rz(-1.1558371) q[3];
sx q[3];
rz(-2.6716148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0726274) q[2];
sx q[2];
rz(-1.5921389) q[2];
sx q[2];
rz(2.3907982) q[2];
rz(0.99610656) q[3];
sx q[3];
rz(-0.80544296) q[3];
sx q[3];
rz(-0.44710818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65425101) q[0];
sx q[0];
rz(-0.60289201) q[0];
sx q[0];
rz(-0.61016369) q[0];
rz(0.24972406) q[1];
sx q[1];
rz(-1.4796673) q[1];
sx q[1];
rz(-3.0009559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2785169) q[0];
sx q[0];
rz(-0.91916537) q[0];
sx q[0];
rz(-0.58653513) q[0];
x q[1];
rz(1.6737559) q[2];
sx q[2];
rz(-2.0741346) q[2];
sx q[2];
rz(-1.9416941) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36547713) q[1];
sx q[1];
rz(-1.2942614) q[1];
sx q[1];
rz(1.4955669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3811627) q[3];
sx q[3];
rz(-1.7685585) q[3];
sx q[3];
rz(-2.9341079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3934624) q[2];
sx q[2];
rz(-2.738214) q[2];
sx q[2];
rz(1.0900137) q[2];
rz(1.9213093) q[3];
sx q[3];
rz(-2.0956495) q[3];
sx q[3];
rz(2.2530341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1423993) q[0];
sx q[0];
rz(-2.0274473) q[0];
sx q[0];
rz(0.92495579) q[0];
rz(1.6389182) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(-2.9138873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3229436) q[0];
sx q[0];
rz(-0.54185757) q[0];
sx q[0];
rz(1.3635554) q[0];
rz(-pi) q[1];
rz(0.57982071) q[2];
sx q[2];
rz(-1.3353773) q[2];
sx q[2];
rz(-1.5583684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0947837) q[1];
sx q[1];
rz(-1.9247479) q[1];
sx q[1];
rz(1.3665529) q[1];
rz(-pi) q[2];
rz(1.3761767) q[3];
sx q[3];
rz(-2.1756919) q[3];
sx q[3];
rz(-2.5218487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0324675) q[2];
sx q[2];
rz(-2.511907) q[2];
sx q[2];
rz(0.49051782) q[2];
rz(2.0971175) q[3];
sx q[3];
rz(-2.677768) q[3];
sx q[3];
rz(-3.131955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57546416) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(1.2658966) q[0];
rz(1.6195541) q[1];
sx q[1];
rz(-1.252389) q[1];
sx q[1];
rz(0.44580805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193008) q[0];
sx q[0];
rz(-2.5133488) q[0];
sx q[0];
rz(-0.99530812) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51165255) q[2];
sx q[2];
rz(-1.5298136) q[2];
sx q[2];
rz(0.042854007) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3109772) q[1];
sx q[1];
rz(-1.4339245) q[1];
sx q[1];
rz(-0.61424436) q[1];
rz(-pi) q[2];
rz(0.021546797) q[3];
sx q[3];
rz(-1.8441219) q[3];
sx q[3];
rz(-0.97238982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.012099115) q[2];
sx q[2];
rz(-0.98126498) q[2];
sx q[2];
rz(-2.1141466) q[2];
rz(1.8217314) q[3];
sx q[3];
rz(-0.3052932) q[3];
sx q[3];
rz(-0.52904883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32578596) q[0];
sx q[0];
rz(-1.3675714) q[0];
sx q[0];
rz(2.0969781) q[0];
rz(-0.927399) q[1];
sx q[1];
rz(-1.8721885) q[1];
sx q[1];
rz(2.4354557) q[1];
rz(1.883222) q[2];
sx q[2];
rz(-1.8835405) q[2];
sx q[2];
rz(0.082688511) q[2];
rz(0.61946034) q[3];
sx q[3];
rz(-1.3347698) q[3];
sx q[3];
rz(2.6556849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
