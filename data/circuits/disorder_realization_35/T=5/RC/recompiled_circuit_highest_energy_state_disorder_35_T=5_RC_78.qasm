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
rz(0.063989446) q[0];
sx q[0];
rz(-2.2292697) q[0];
sx q[0];
rz(1.8066701) q[0];
rz(1.9995243) q[1];
sx q[1];
rz(-2.6231782) q[1];
sx q[1];
rz(-1.6496744) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.175615) q[0];
sx q[0];
rz(-1.5892995) q[0];
sx q[0];
rz(-1.6241839) q[0];
x q[1];
rz(0.71662997) q[2];
sx q[2];
rz(-2.306005) q[2];
sx q[2];
rz(0.40975964) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8644997) q[1];
sx q[1];
rz(-0.96597176) q[1];
sx q[1];
rz(0.34326174) q[1];
x q[2];
rz(-1.8133477) q[3];
sx q[3];
rz(-1.6036766) q[3];
sx q[3];
rz(-1.3803079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(-0.054923687) q[2];
rz(-1.6867636) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(-1.8964881) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354483) q[0];
sx q[0];
rz(-1.7559253) q[0];
sx q[0];
rz(-1.5371176) q[0];
x q[1];
rz(-2.6284559) q[2];
sx q[2];
rz(-1.0982795) q[2];
sx q[2];
rz(-3.0233011) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5692484) q[1];
sx q[1];
rz(-0.73671417) q[1];
sx q[1];
rz(2.8235497) q[1];
rz(-pi) q[2];
rz(3.008083) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(0.51160073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9241141) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(-0.58383101) q[2];
rz(-0.42932388) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(-2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15447021) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(-1.467147) q[0];
rz(-1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(1.0999701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3008793) q[0];
sx q[0];
rz(-1.3726808) q[0];
sx q[0];
rz(-1.677657) q[0];
rz(-pi) q[1];
rz(1.9103545) q[2];
sx q[2];
rz(-1.8446088) q[2];
sx q[2];
rz(-1.6664315) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63640672) q[1];
sx q[1];
rz(-0.51915324) q[1];
sx q[1];
rz(0.7208419) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6925642) q[3];
sx q[3];
rz(-2.0981023) q[3];
sx q[3];
rz(-2.9849986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35525068) q[2];
sx q[2];
rz(-0.92266551) q[2];
sx q[2];
rz(1.5477017) q[2];
rz(1.8111546) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(-1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610157) q[0];
sx q[0];
rz(-1.3207734) q[0];
sx q[0];
rz(2.9837578) q[0];
rz(0.24179587) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(0.48113021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24836981) q[0];
sx q[0];
rz(-0.71740323) q[0];
sx q[0];
rz(-0.1347085) q[0];
rz(2.782576) q[2];
sx q[2];
rz(-1.8098157) q[2];
sx q[2];
rz(-1.1306819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9485875) q[1];
sx q[1];
rz(-0.95110106) q[1];
sx q[1];
rz(-2.7189458) q[1];
rz(-pi) q[2];
rz(2.1948482) q[3];
sx q[3];
rz(-2.7928154) q[3];
sx q[3];
rz(2.6117976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87159291) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(0.2743741) q[2];
rz(-1.4659878) q[3];
sx q[3];
rz(-1.8612739) q[3];
sx q[3];
rz(2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(-1.0614606) q[0];
rz(0.44867107) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(-1.7015069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7913623) q[0];
sx q[0];
rz(-1.7273979) q[0];
sx q[0];
rz(-0.76275218) q[0];
rz(-3.0147047) q[2];
sx q[2];
rz(-1.1694307) q[2];
sx q[2];
rz(0.021406476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49855194) q[1];
sx q[1];
rz(-1.3128377) q[1];
sx q[1];
rz(3.113913) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1813356) q[3];
sx q[3];
rz(-1.6597851) q[3];
sx q[3];
rz(-1.2582092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3638641) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(-2.7423972) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6804009) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(0.68049085) q[0];
rz(2.3091799) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(2.9852273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6844821) q[0];
sx q[0];
rz(-1.5183987) q[0];
sx q[0];
rz(0.27586097) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83145468) q[2];
sx q[2];
rz(-2.105684) q[2];
sx q[2];
rz(-1.8597459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.029441) q[1];
sx q[1];
rz(-1.8990771) q[1];
sx q[1];
rz(2.8480457) q[1];
x q[2];
rz(-0.12229192) q[3];
sx q[3];
rz(-1.3548791) q[3];
sx q[3];
rz(-1.9401039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6303595) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(0.88999256) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(0.61221468) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(0.30119687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071547) q[0];
sx q[0];
rz(-0.2858735) q[0];
sx q[0];
rz(1.8854333) q[0];
x q[1];
rz(-2.7686504) q[2];
sx q[2];
rz(-0.94764793) q[2];
sx q[2];
rz(2.9187893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9928678) q[1];
sx q[1];
rz(-0.5827671) q[1];
sx q[1];
rz(0.73281835) q[1];
rz(1.3314358) q[3];
sx q[3];
rz(-1.3151957) q[3];
sx q[3];
rz(-1.1472697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(2.0937199) q[2];
rz(2.786484) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(-2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316786) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-0.32383305) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-2.8388265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51311868) q[0];
sx q[0];
rz(-1.895005) q[0];
sx q[0];
rz(-1.8379492) q[0];
rz(-pi) q[1];
rz(-2.9043973) q[2];
sx q[2];
rz(-1.5233808) q[2];
sx q[2];
rz(-0.42979017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0787133) q[1];
sx q[1];
rz(-2.4643366) q[1];
sx q[1];
rz(0.96911624) q[1];
rz(-pi) q[2];
rz(1.15074) q[3];
sx q[3];
rz(-1.8057293) q[3];
sx q[3];
rz(0.7701503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6451463) q[2];
sx q[2];
rz(-1.030693) q[2];
sx q[2];
rz(0.56078625) q[2];
rz(2.136611) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(-1.0952449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.1110558) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(-0.38052446) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(-2.1662625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935287) q[0];
sx q[0];
rz(-0.87665999) q[0];
sx q[0];
rz(2.3175236) q[0];
rz(-pi) q[1];
rz(-3.0815711) q[2];
sx q[2];
rz(-1.2916471) q[2];
sx q[2];
rz(0.44250968) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0142648) q[1];
sx q[1];
rz(-1.2593966) q[1];
sx q[1];
rz(2.6173956) q[1];
rz(-pi) q[2];
rz(0.22244979) q[3];
sx q[3];
rz(-1.8651903) q[3];
sx q[3];
rz(-2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9425977) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(-1.0814166) q[2];
rz(-0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-2.2190905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0104495) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(-2.8358054) q[0];
rz(2.8336613) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(1.9889132) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891805) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(-2.7680725) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20736097) q[2];
sx q[2];
rz(-2.2054858) q[2];
sx q[2];
rz(0.35516741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.071142405) q[1];
sx q[1];
rz(-1.2463039) q[1];
sx q[1];
rz(-1.9981415) q[1];
rz(-1.6341428) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(1.9346332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(1.6884241) q[2];
rz(2.6609663) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(-1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4881445) q[0];
sx q[0];
rz(-1.6076037) q[0];
sx q[0];
rz(-1.6758767) q[0];
rz(-1.0427955) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(0.53521143) q[2];
sx q[2];
rz(-1.3387398) q[2];
sx q[2];
rz(-2.6571318) q[2];
rz(-1.5299464) q[3];
sx q[3];
rz(-1.9286641) q[3];
sx q[3];
rz(2.7947938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
