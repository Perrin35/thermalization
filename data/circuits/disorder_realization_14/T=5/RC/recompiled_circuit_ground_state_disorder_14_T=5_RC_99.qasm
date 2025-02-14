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
rz(-0.39363632) q[1];
sx q[1];
rz(4.7596158) q[1];
sx q[1];
rz(14.16852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868139) q[0];
sx q[0];
rz(-0.20113763) q[0];
sx q[0];
rz(3.1166409) q[0];
rz(0.27199409) q[2];
sx q[2];
rz(-1.0390721) q[2];
sx q[2];
rz(3.0716999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56583079) q[1];
sx q[1];
rz(-1.1850433) q[1];
sx q[1];
rz(-1.7931531) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39253916) q[3];
sx q[3];
rz(-1.2219567) q[3];
sx q[3];
rz(2.4667645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9615122) q[2];
sx q[2];
rz(-0.95246685) q[2];
sx q[2];
rz(-0.81641475) q[2];
rz(-0.69711971) q[3];
sx q[3];
rz(-0.37728798) q[3];
sx q[3];
rz(-0.94059801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78696752) q[0];
sx q[0];
rz(-2.0836199) q[0];
sx q[0];
rz(0.70607287) q[0];
rz(-0.1708897) q[1];
sx q[1];
rz(-0.51097521) q[1];
sx q[1];
rz(-0.99542803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1814897) q[0];
sx q[0];
rz(-1.2523635) q[0];
sx q[0];
rz(-2.7820827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.273492) q[2];
sx q[2];
rz(-1.9579525) q[2];
sx q[2];
rz(2.0086945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1198955) q[1];
sx q[1];
rz(-1.3576389) q[1];
sx q[1];
rz(-2.9851856) q[1];
rz(-pi) q[2];
rz(-1.6042406) q[3];
sx q[3];
rz(-2.20801) q[3];
sx q[3];
rz(2.6496365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.55676111) q[2];
sx q[2];
rz(-2.2274667) q[2];
sx q[2];
rz(-2.8625028) q[2];
rz(2.4537405) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85182178) q[0];
sx q[0];
rz(-0.91854799) q[0];
sx q[0];
rz(-0.37891349) q[0];
rz(-1.9508349) q[1];
sx q[1];
rz(-2.7749116) q[1];
sx q[1];
rz(0.77622882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8541753) q[0];
sx q[0];
rz(-2.6834134) q[0];
sx q[0];
rz(-1.2121546) q[0];
rz(-pi) q[1];
rz(0.68910901) q[2];
sx q[2];
rz(-2.9053134) q[2];
sx q[2];
rz(-0.29321996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32404985) q[1];
sx q[1];
rz(-0.88081283) q[1];
sx q[1];
rz(1.2381366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.046136304) q[3];
sx q[3];
rz(-2.0946119) q[3];
sx q[3];
rz(0.55098096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0188521) q[2];
sx q[2];
rz(-1.2383702) q[2];
sx q[2];
rz(-2.5416601) q[2];
rz(-1.9020724) q[3];
sx q[3];
rz(-1.4116838) q[3];
sx q[3];
rz(-0.6111353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3032853) q[0];
sx q[0];
rz(-0.67876434) q[0];
sx q[0];
rz(-1.7093866) q[0];
rz(0.86668658) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(-2.0492882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0129654) q[0];
sx q[0];
rz(-2.1218637) q[0];
sx q[0];
rz(2.4222213) q[0];
rz(-pi) q[1];
rz(2.1020911) q[2];
sx q[2];
rz(-0.73248011) q[2];
sx q[2];
rz(1.3074444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6505423) q[1];
sx q[1];
rz(-2.67727) q[1];
sx q[1];
rz(3.1373255) q[1];
rz(2.3471429) q[3];
sx q[3];
rz(-0.91728079) q[3];
sx q[3];
rz(0.78032035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0852574) q[2];
sx q[2];
rz(-1.9838355) q[2];
sx q[2];
rz(-1.2752656) q[2];
rz(1.3750252) q[3];
sx q[3];
rz(-1.2396038) q[3];
sx q[3];
rz(1.0331155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.105724) q[0];
sx q[0];
rz(-0.88304702) q[0];
sx q[0];
rz(1.7778273) q[0];
rz(-0.56963244) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(2.6380576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1138332) q[0];
sx q[0];
rz(-0.22969023) q[0];
sx q[0];
rz(-2.5656347) q[0];
x q[1];
rz(0.72959186) q[2];
sx q[2];
rz(-0.9895095) q[2];
sx q[2];
rz(-0.94129291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9787748) q[1];
sx q[1];
rz(-1.4895338) q[1];
sx q[1];
rz(0.012465076) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3713417) q[3];
sx q[3];
rz(-1.7903084) q[3];
sx q[3];
rz(1.3971922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2786431) q[2];
sx q[2];
rz(-2.7549665) q[2];
sx q[2];
rz(1.53842) q[2];
rz(-3.0053511) q[3];
sx q[3];
rz(-1.7377661) q[3];
sx q[3];
rz(-1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513026) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(2.8001617) q[0];
rz(2.4104207) q[1];
sx q[1];
rz(-2.5656504) q[1];
sx q[1];
rz(2.0700571) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2444457) q[0];
sx q[0];
rz(-1.1700774) q[0];
sx q[0];
rz(2.6942376) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5265556) q[2];
sx q[2];
rz(-1.5205548) q[2];
sx q[2];
rz(3.0453403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4298297) q[1];
sx q[1];
rz(-2.2380658) q[1];
sx q[1];
rz(-2.1416452) q[1];
rz(-pi) q[2];
rz(0.0045147379) q[3];
sx q[3];
rz(-1.0663973) q[3];
sx q[3];
rz(0.90108353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3557446) q[2];
sx q[2];
rz(-1.5037437) q[2];
sx q[2];
rz(2.6453633) q[2];
rz(-1.5931386) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(-2.1350071) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9736495) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(3.1267401) q[0];
rz(1.4133833) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(-2.3766439) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1122862) q[0];
sx q[0];
rz(-2.1878982) q[0];
sx q[0];
rz(-0.86227148) q[0];
x q[1];
rz(1.9221481) q[2];
sx q[2];
rz(-0.91809154) q[2];
sx q[2];
rz(1.2437133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98337631) q[1];
sx q[1];
rz(-1.9675273) q[1];
sx q[1];
rz(-1.5660172) q[1];
rz(0.21723218) q[3];
sx q[3];
rz(-1.1558371) q[3];
sx q[3];
rz(2.6716148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0689653) q[2];
sx q[2];
rz(-1.5921389) q[2];
sx q[2];
rz(0.7507945) q[2];
rz(-0.99610656) q[3];
sx q[3];
rz(-0.80544296) q[3];
sx q[3];
rz(-2.6944845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.4873416) q[0];
sx q[0];
rz(-0.60289201) q[0];
sx q[0];
rz(0.61016369) q[0];
rz(0.24972406) q[1];
sx q[1];
rz(-1.6619253) q[1];
sx q[1];
rz(-0.14063674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8630757) q[0];
sx q[0];
rz(-0.91916537) q[0];
sx q[0];
rz(2.5550575) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9570661) q[2];
sx q[2];
rz(-2.6287196) q[2];
sx q[2];
rz(-1.41092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0454694) q[1];
sx q[1];
rz(-2.8552607) q[1];
sx q[1];
rz(-0.25888049) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.300977) q[3];
sx q[3];
rz(-0.82874819) q[3];
sx q[3];
rz(-1.5935115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3934624) q[2];
sx q[2];
rz(-2.738214) q[2];
sx q[2];
rz(-2.051579) q[2];
rz(1.2202834) q[3];
sx q[3];
rz(-2.0956495) q[3];
sx q[3];
rz(0.88855851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1423993) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(-2.2166369) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(-0.22770539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.071969) q[0];
sx q[0];
rz(-1.6771131) q[0];
sx q[0];
rz(-2.1031455) q[0];
rz(-0.57982071) q[2];
sx q[2];
rz(-1.3353773) q[2];
sx q[2];
rz(-1.5832242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.633435) q[1];
sx q[1];
rz(-2.7351008) q[1];
sx q[1];
rz(0.50199957) q[1];
x q[2];
rz(1.7654159) q[3];
sx q[3];
rz(-0.96590079) q[3];
sx q[3];
rz(-2.5218487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0324675) q[2];
sx q[2];
rz(-2.511907) q[2];
sx q[2];
rz(-2.6510748) q[2];
rz(-1.0444752) q[3];
sx q[3];
rz(-2.677768) q[3];
sx q[3];
rz(0.0096376816) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5661285) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(1.2658966) q[0];
rz(1.6195541) q[1];
sx q[1];
rz(-1.252389) q[1];
sx q[1];
rz(-2.6957846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193008) q[0];
sx q[0];
rz(-2.5133488) q[0];
sx q[0];
rz(-0.99530812) q[0];
rz(-pi) q[1];
rz(-1.61779) q[2];
sx q[2];
rz(-2.0819774) q[2];
sx q[2];
rz(1.5509449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3054786) q[1];
sx q[1];
rz(-2.1784601) q[1];
sx q[1];
rz(1.4038248) q[1];
rz(-pi) q[2];
x q[2];
rz(0.021546797) q[3];
sx q[3];
rz(-1.2974707) q[3];
sx q[3];
rz(0.97238982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.012099115) q[2];
sx q[2];
rz(-0.98126498) q[2];
sx q[2];
rz(2.1141466) q[2];
rz(1.3198613) q[3];
sx q[3];
rz(-2.8362995) q[3];
sx q[3];
rz(-0.52904883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8158067) q[0];
sx q[0];
rz(-1.7740213) q[0];
sx q[0];
rz(-1.0446145) q[0];
rz(2.2141937) q[1];
sx q[1];
rz(-1.8721885) q[1];
sx q[1];
rz(2.4354557) q[1];
rz(-0.32756424) q[2];
sx q[2];
rz(-1.8675928) q[2];
sx q[2];
rz(1.7525385) q[2];
rz(0.3927265) q[3];
sx q[3];
rz(-0.65734335) q[3];
sx q[3];
rz(-1.7398294) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
