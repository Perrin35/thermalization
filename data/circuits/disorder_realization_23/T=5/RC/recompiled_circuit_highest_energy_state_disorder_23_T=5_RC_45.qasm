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
rz(0.55573207) q[0];
sx q[0];
rz(-1.8620123) q[0];
sx q[0];
rz(-0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965138) q[0];
sx q[0];
rz(-0.9849087) q[0];
sx q[0];
rz(-1.3784598) q[0];
rz(-1.3921803) q[2];
sx q[2];
rz(-1.1682604) q[2];
sx q[2];
rz(0.91031633) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7465621) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(-1.0814352) q[1];
x q[2];
rz(-1.1527083) q[3];
sx q[3];
rz(-1.3527414) q[3];
sx q[3];
rz(-2.5421028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27935394) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(-1.5581101) q[2];
rz(-2.8067348) q[3];
sx q[3];
rz(-1.0840651) q[3];
sx q[3];
rz(-2.8607232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25664499) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(0.33194342) q[0];
rz(2.7815107) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(-2.8538381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641787) q[0];
sx q[0];
rz(-2.8197643) q[0];
sx q[0];
rz(-0.68049707) q[0];
rz(-pi) q[1];
rz(0.73189484) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(-0.11920028) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1102241) q[1];
sx q[1];
rz(-1.2249631) q[1];
sx q[1];
rz(-0.35067534) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5310181) q[3];
sx q[3];
rz(-2.861851) q[3];
sx q[3];
rz(0.66204643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45431367) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(1.3909371) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39621064) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(-2.9111653) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968937) q[0];
sx q[0];
rz(-2.7285577) q[0];
sx q[0];
rz(1.4223736) q[0];
rz(-pi) q[1];
rz(2.6730113) q[2];
sx q[2];
rz(-1.5043253) q[2];
sx q[2];
rz(-2.2867212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9625712) q[1];
sx q[1];
rz(-1.6238535) q[1];
sx q[1];
rz(-1.9808484) q[1];
rz(-2.2346441) q[3];
sx q[3];
rz(-2.8606599) q[3];
sx q[3];
rz(-0.30981608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0383703) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(3.1008516) q[2];
rz(3.0401958) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(1.066677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7539702) q[0];
sx q[0];
rz(-0.0059703537) q[0];
sx q[0];
rz(2.907584) q[0];
rz(-0.19800828) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(-2.1580946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6925544) q[0];
sx q[0];
rz(-0.48547599) q[0];
sx q[0];
rz(-1.4794769) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6621236) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(0.78081607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9069104) q[1];
sx q[1];
rz(-2.249472) q[1];
sx q[1];
rz(-0.80295103) q[1];
x q[2];
rz(-1.6487021) q[3];
sx q[3];
rz(-2.7512105) q[3];
sx q[3];
rz(1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46821758) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(-0.77486983) q[2];
rz(-1.3261999) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(2.6302122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.4020017) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(0.032489754) q[0];
rz(1.912311) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(0.083018735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0262526) q[0];
sx q[0];
rz(-1.3466382) q[0];
sx q[0];
rz(0.46333939) q[0];
rz(2.584721) q[2];
sx q[2];
rz(-0.31302127) q[2];
sx q[2];
rz(1.0525296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0159576) q[1];
sx q[1];
rz(-2.2803366) q[1];
sx q[1];
rz(-0.63035359) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.015542726) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(1.6890845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7897537) q[2];
sx q[2];
rz(-0.80078501) q[2];
sx q[2];
rz(-2.9808673) q[2];
rz(-1.9488581) q[3];
sx q[3];
rz(-2.9294117) q[3];
sx q[3];
rz(2.3737657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47650325) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(-0.29092586) q[0];
rz(1.3833969) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(2.6417522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0394589) q[0];
sx q[0];
rz(-1.5284662) q[0];
sx q[0];
rz(-3.0979587) q[0];
rz(-pi) q[1];
rz(-0.84264522) q[2];
sx q[2];
rz(-2.004874) q[2];
sx q[2];
rz(1.7378716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0840869) q[1];
sx q[1];
rz(-2.018714) q[1];
sx q[1];
rz(0.65051669) q[1];
rz(-pi) q[2];
rz(-0.96792241) q[3];
sx q[3];
rz(-0.80652666) q[3];
sx q[3];
rz(1.4462808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2073888) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(-0.70154166) q[2];
rz(-0.97405854) q[3];
sx q[3];
rz(-0.8166703) q[3];
sx q[3];
rz(-2.2402703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.3153673) q[0];
sx q[0];
rz(-2.3235445) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(1.2391799) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(2.0674131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095431134) q[0];
sx q[0];
rz(-0.019733075) q[0];
sx q[0];
rz(-2.31953) q[0];
x q[1];
rz(-2.7137382) q[2];
sx q[2];
rz(-1.9249467) q[2];
sx q[2];
rz(-0.71989518) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3201808) q[1];
sx q[1];
rz(-1.786725) q[1];
sx q[1];
rz(1.126299) q[1];
x q[2];
rz(-2.6857008) q[3];
sx q[3];
rz(-2.3237202) q[3];
sx q[3];
rz(-2.4867833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2913975) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(-0.58471739) q[2];
rz(-2.0753453) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(-0.48072746) q[0];
rz(1.9302543) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(-0.617625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4942426) q[0];
sx q[0];
rz(-1.8959008) q[0];
sx q[0];
rz(-0.0044429739) q[0];
rz(-0.60725339) q[2];
sx q[2];
rz(-1.1285845) q[2];
sx q[2];
rz(0.84458015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94719515) q[1];
sx q[1];
rz(-0.36809599) q[1];
sx q[1];
rz(-0.93666623) q[1];
rz(1.6061574) q[3];
sx q[3];
rz(-1.3786982) q[3];
sx q[3];
rz(1.1079605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(0.89087957) q[2];
rz(1.0293055) q[3];
sx q[3];
rz(-1.7512713) q[3];
sx q[3];
rz(-1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35128281) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(-1.9679605) q[0];
rz(-1.1890746) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(0.48114166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3404589) q[0];
sx q[0];
rz(-2.1251571) q[0];
sx q[0];
rz(1.3075243) q[0];
x q[1];
rz(1.6716154) q[2];
sx q[2];
rz(-2.1283796) q[2];
sx q[2];
rz(2.4876311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18080841) q[1];
sx q[1];
rz(-0.97429064) q[1];
sx q[1];
rz(2.9636613) q[1];
rz(-pi) q[2];
rz(-1.0403544) q[3];
sx q[3];
rz(-1.9602437) q[3];
sx q[3];
rz(-2.3688909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(-1.6602328) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.1888844) q[3];
sx q[3];
rz(-1.9143298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34761053) q[0];
sx q[0];
rz(-2.1078258) q[0];
sx q[0];
rz(0.069742918) q[0];
rz(0.28930411) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(-2.0893673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7523868) q[0];
sx q[0];
rz(-1.9294039) q[0];
sx q[0];
rz(-1.9176656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90862008) q[2];
sx q[2];
rz(-1.7025456) q[2];
sx q[2];
rz(0.8140623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94967647) q[1];
sx q[1];
rz(-1.9468309) q[1];
sx q[1];
rz(-2.7533349) q[1];
x q[2];
rz(-1.9931562) q[3];
sx q[3];
rz(-1.1760654) q[3];
sx q[3];
rz(2.0572061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5436486) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(-0.0607461) q[2];
rz(2.8525823) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(-1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48988265) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(-1.0501077) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(-2.9906828) q[2];
sx q[2];
rz(-1.1058634) q[2];
sx q[2];
rz(2.4185218) q[2];
rz(-1.0084739) q[3];
sx q[3];
rz(-0.49839603) q[3];
sx q[3];
rz(-2.8965542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
