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
rz(0.93208575) q[0];
sx q[0];
rz(-2.1075489) q[0];
sx q[0];
rz(2.3873868) q[0];
rz(1.7475313) q[1];
sx q[1];
rz(-0.59352195) q[1];
sx q[1];
rz(-1.4374179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0707929) q[0];
sx q[0];
rz(-1.6294998) q[0];
sx q[0];
rz(-3.015201) q[0];
rz(-pi) q[1];
rz(-1.9600771) q[2];
sx q[2];
rz(-1.5354275) q[2];
sx q[2];
rz(-0.25098342) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5240548) q[1];
sx q[1];
rz(-1.0144985) q[1];
sx q[1];
rz(1.9281045) q[1];
rz(-pi) q[2];
rz(3.0079477) q[3];
sx q[3];
rz(-1.9507532) q[3];
sx q[3];
rz(-0.28912543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6393911) q[2];
sx q[2];
rz(-0.23362564) q[2];
sx q[2];
rz(0.26546738) q[2];
rz(1.4207077) q[3];
sx q[3];
rz(-1.1666433) q[3];
sx q[3];
rz(1.7336806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.0625286) q[0];
sx q[0];
rz(-1.9363576) q[0];
sx q[0];
rz(-1.0519387) q[0];
rz(0.99986347) q[1];
sx q[1];
rz(-1.544416) q[1];
sx q[1];
rz(2.3725407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3628415) q[0];
sx q[0];
rz(-1.0386416) q[0];
sx q[0];
rz(0.51409419) q[0];
rz(-pi) q[1];
rz(1.6300419) q[2];
sx q[2];
rz(-0.99397221) q[2];
sx q[2];
rz(1.8920997) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1726672) q[1];
sx q[1];
rz(-1.7255867) q[1];
sx q[1];
rz(-1.118578) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7160109) q[3];
sx q[3];
rz(-1.3481529) q[3];
sx q[3];
rz(-2.1964335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96180463) q[2];
sx q[2];
rz(-2.3369117) q[2];
sx q[2];
rz(1.556832) q[2];
rz(2.4368317) q[3];
sx q[3];
rz(-0.2125936) q[3];
sx q[3];
rz(1.0889277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8162808) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(2.3749206) q[0];
rz(2.7491772) q[1];
sx q[1];
rz(-2.2800443) q[1];
sx q[1];
rz(0.79889417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2080059) q[0];
sx q[0];
rz(-2.0547935) q[0];
sx q[0];
rz(0.010864929) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8908752) q[2];
sx q[2];
rz(-1.6581767) q[2];
sx q[2];
rz(2.6911497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0392363) q[1];
sx q[1];
rz(-1.2627541) q[1];
sx q[1];
rz(2.0935489) q[1];
rz(0.17744448) q[3];
sx q[3];
rz(-1.6758306) q[3];
sx q[3];
rz(0.80897409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40156349) q[2];
sx q[2];
rz(-2.8758958) q[2];
sx q[2];
rz(-3.1166039) q[2];
rz(-1.7449024) q[3];
sx q[3];
rz(-1.2324421) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5806737) q[0];
sx q[0];
rz(-2.6011401) q[0];
sx q[0];
rz(0.29676357) q[0];
rz(0.98776039) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(0.74772778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4241705) q[0];
sx q[0];
rz(-1.7892196) q[0];
sx q[0];
rz(2.3481247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.426366) q[2];
sx q[2];
rz(-1.919121) q[2];
sx q[2];
rz(2.0722318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8104483) q[1];
sx q[1];
rz(-0.40256009) q[1];
sx q[1];
rz(-0.029235201) q[1];
x q[2];
rz(0.70246299) q[3];
sx q[3];
rz(-0.80804658) q[3];
sx q[3];
rz(-1.3987227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1838386) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(2.6217065) q[2];
rz(-0.94684354) q[3];
sx q[3];
rz(-1.9498473) q[3];
sx q[3];
rz(0.051518353) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87676048) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-1.1967891) q[0];
rz(1.4747249) q[1];
sx q[1];
rz(-0.80497545) q[1];
sx q[1];
rz(2.8428452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099649053) q[0];
sx q[0];
rz(-1.9199326) q[0];
sx q[0];
rz(0.70167244) q[0];
rz(2.0709956) q[2];
sx q[2];
rz(-0.83865863) q[2];
sx q[2];
rz(-0.054890779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1890002) q[1];
sx q[1];
rz(-2.3680284) q[1];
sx q[1];
rz(2.9783334) q[1];
x q[2];
rz(2.7389333) q[3];
sx q[3];
rz(-2.23051) q[3];
sx q[3];
rz(-0.92056489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8220736) q[2];
sx q[2];
rz(-0.75054979) q[2];
sx q[2];
rz(0.015241148) q[2];
rz(-3.0139253) q[3];
sx q[3];
rz(-2.7805507) q[3];
sx q[3];
rz(2.291919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.78758883) q[0];
sx q[0];
rz(-0.51743999) q[0];
sx q[0];
rz(2.2678243) q[0];
rz(-1.9173701) q[1];
sx q[1];
rz(-1.0226378) q[1];
sx q[1];
rz(1.9220985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2827493) q[0];
sx q[0];
rz(-1.97012) q[0];
sx q[0];
rz(2.0612677) q[0];
rz(-0.60004931) q[2];
sx q[2];
rz(-0.82227899) q[2];
sx q[2];
rz(-2.6278969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.022910206) q[1];
sx q[1];
rz(-1.3617651) q[1];
sx q[1];
rz(0.19281705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9683011) q[3];
sx q[3];
rz(-2.1263206) q[3];
sx q[3];
rz(1.3501957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6957317) q[2];
sx q[2];
rz(-1.4706688) q[2];
sx q[2];
rz(-0.50103465) q[2];
rz(-1.3387574) q[3];
sx q[3];
rz(-2.7300291) q[3];
sx q[3];
rz(1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434791) q[0];
sx q[0];
rz(-1.1971594) q[0];
sx q[0];
rz(-0.57394779) q[0];
rz(-0.57169882) q[1];
sx q[1];
rz(-1.0400925) q[1];
sx q[1];
rz(1.5843102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27423647) q[0];
sx q[0];
rz(-0.78834962) q[0];
sx q[0];
rz(-1.0324182) q[0];
x q[1];
rz(-0.90960501) q[2];
sx q[2];
rz(-1.8665458) q[2];
sx q[2];
rz(-0.32509229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3347335) q[1];
sx q[1];
rz(-0.99913952) q[1];
sx q[1];
rz(0.57493361) q[1];
rz(0.81767003) q[3];
sx q[3];
rz(-0.51556057) q[3];
sx q[3];
rz(1.8918693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68465662) q[2];
sx q[2];
rz(-1.9921649) q[2];
sx q[2];
rz(-1.9880902) q[2];
rz(1.614511) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1077147) q[0];
sx q[0];
rz(-2.3623473) q[0];
sx q[0];
rz(0.23767924) q[0];
rz(2.4029845) q[1];
sx q[1];
rz(-1.9269201) q[1];
sx q[1];
rz(-0.019066378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8066609) q[0];
sx q[0];
rz(-1.4981759) q[0];
sx q[0];
rz(-0.0082834335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7007199) q[2];
sx q[2];
rz(-2.2402175) q[2];
sx q[2];
rz(0.6159516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70266378) q[1];
sx q[1];
rz(-3.0512718) q[1];
sx q[1];
rz(2.0144736) q[1];
rz(-pi) q[2];
rz(0.98293368) q[3];
sx q[3];
rz(-1.0140099) q[3];
sx q[3];
rz(0.24036053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1261403) q[2];
sx q[2];
rz(-1.3498053) q[2];
sx q[2];
rz(0.17507412) q[2];
rz(1.3779047) q[3];
sx q[3];
rz(-2.521069) q[3];
sx q[3];
rz(0.11848816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058067583) q[0];
sx q[0];
rz(-1.7836934) q[0];
sx q[0];
rz(2.626626) q[0];
rz(-0.99634755) q[1];
sx q[1];
rz(-1.0641655) q[1];
sx q[1];
rz(1.4774342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2124951) q[0];
sx q[0];
rz(-0.97508865) q[0];
sx q[0];
rz(-0.86796772) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58837946) q[2];
sx q[2];
rz(-1.5984238) q[2];
sx q[2];
rz(-0.80828062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2075726) q[1];
sx q[1];
rz(-2.0193384) q[1];
sx q[1];
rz(2.1357082) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1636073) q[3];
sx q[3];
rz(-1.8043274) q[3];
sx q[3];
rz(1.1790891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.065206334) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(0.15929407) q[2];
rz(-1.6431036) q[3];
sx q[3];
rz(-2.4703333) q[3];
sx q[3];
rz(1.8026277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945187) q[0];
sx q[0];
rz(-0.044059489) q[0];
sx q[0];
rz(-0.86167589) q[0];
rz(1.8623976) q[1];
sx q[1];
rz(-2.8586614) q[1];
sx q[1];
rz(0.18712015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30062723) q[0];
sx q[0];
rz(-1.297313) q[0];
sx q[0];
rz(1.8086368) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2995903) q[2];
sx q[2];
rz(-2.0496297) q[2];
sx q[2];
rz(1.9128654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5622065) q[1];
sx q[1];
rz(-2.36464) q[1];
sx q[1];
rz(0.38511606) q[1];
rz(-pi) q[2];
rz(-2.7469603) q[3];
sx q[3];
rz(-1.8642211) q[3];
sx q[3];
rz(-0.59673264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44237915) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(1.7787735) q[2];
rz(-1.5022701) q[3];
sx q[3];
rz(-0.5539186) q[3];
sx q[3];
rz(1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61158553) q[0];
sx q[0];
rz(-1.3349608) q[0];
sx q[0];
rz(0.0064749574) q[0];
rz(0.50637983) q[1];
sx q[1];
rz(-0.19229278) q[1];
sx q[1];
rz(-2.2813588) q[1];
rz(-0.14534563) q[2];
sx q[2];
rz(-1.9898562) q[2];
sx q[2];
rz(-2.936292) q[2];
rz(-1.3527444) q[3];
sx q[3];
rz(-1.8321804) q[3];
sx q[3];
rz(1.2294168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
