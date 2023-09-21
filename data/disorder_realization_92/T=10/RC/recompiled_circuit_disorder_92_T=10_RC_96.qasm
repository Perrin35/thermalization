OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(3.175088) q[0];
sx q[0];
rz(7.6498084) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(2.0096013) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(2.3495673) q[0];
rz(-pi) q[1];
rz(-1.4049454) q[2];
sx q[2];
rz(-1.7300786) q[2];
sx q[2];
rz(2.7878891) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34239009) q[1];
sx q[1];
rz(-2.4203165) q[1];
sx q[1];
rz(0.93625416) q[1];
x q[2];
rz(1.8559998) q[3];
sx q[3];
rz(-1.0343026) q[3];
sx q[3];
rz(-0.59629089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2312317) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3020878) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(0.50305811) q[0];
rz(1.5548276) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(-0.15393004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0921558) q[0];
sx q[0];
rz(-1.0679809) q[0];
sx q[0];
rz(-0.016331971) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2203127) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(2.9608179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.82984) q[1];
sx q[1];
rz(-2.0318188) q[1];
sx q[1];
rz(2.6677569) q[1];
rz(-2.6318195) q[3];
sx q[3];
rz(-1.5715277) q[3];
sx q[3];
rz(2.7754806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.8187693) q[2];
rz(1.4860738) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948497) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(0.90993607) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(2.1562703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6353778) q[0];
sx q[0];
rz(-2.2197476) q[0];
sx q[0];
rz(-0.20091591) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91395949) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(2.7514806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9039771) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(-1.6275089) q[1];
rz(0.91089532) q[3];
sx q[3];
rz(-1.0064555) q[3];
sx q[3];
rz(-2.6499555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2788006) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(0.4425846) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(0.32143337) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-1.1234269) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.4588979) q[1];
sx q[1];
rz(1.7480063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14560315) q[0];
sx q[0];
rz(-1.9847426) q[0];
sx q[0];
rz(-0.67514174) q[0];
rz(-1.7384999) q[2];
sx q[2];
rz(-2.0003013) q[2];
sx q[2];
rz(1.475856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2323944) q[1];
sx q[1];
rz(-2.0226139) q[1];
sx q[1];
rz(2.3934445) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23354236) q[3];
sx q[3];
rz(-0.84905784) q[3];
sx q[3];
rz(-0.28802179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0397772) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(-0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(-2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(-1.2258688) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(-1.3668758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0492832) q[0];
sx q[0];
rz(-1.7607848) q[0];
sx q[0];
rz(2.7490194) q[0];
x q[1];
rz(1.5405802) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(-3.1291762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1495684) q[1];
sx q[1];
rz(-0.7115041) q[1];
sx q[1];
rz(-1.6929388) q[1];
rz(1.9277875) q[3];
sx q[3];
rz(-1.095466) q[3];
sx q[3];
rz(0.19015042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(0.40536353) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6500403) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(0.50338411) q[0];
rz(2.9227496) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(2.4898081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359905) q[0];
sx q[0];
rz(-0.30797568) q[0];
sx q[0];
rz(0.67291798) q[0];
x q[1];
rz(1.4621554) q[2];
sx q[2];
rz(-0.75192736) q[2];
sx q[2];
rz(-1.3736563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3411322) q[1];
sx q[1];
rz(-0.144185) q[1];
sx q[1];
rz(-1.202281) q[1];
rz(-2.7995293) q[3];
sx q[3];
rz(-0.29933375) q[3];
sx q[3];
rz(2.743486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(2.9351249) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(2.5580653) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(-3.0113509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1657432) q[0];
sx q[0];
rz(-1.2184869) q[0];
sx q[0];
rz(1.6238814) q[0];
rz(-pi) q[1];
rz(0.98353705) q[2];
sx q[2];
rz(-0.98896719) q[2];
sx q[2];
rz(0.26435095) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0252467) q[1];
sx q[1];
rz(-1.3533918) q[1];
sx q[1];
rz(-2.6245963) q[1];
rz(2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(0.91528085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6234201) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(2.0557892) q[2];
rz(-0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(2.3988147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966492) q[0];
sx q[0];
rz(-1.3784694) q[0];
sx q[0];
rz(-2.3924475) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0413264) q[2];
sx q[2];
rz(-0.96457446) q[2];
sx q[2];
rz(0.60339061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10272664) q[1];
sx q[1];
rz(-2.3042149) q[1];
sx q[1];
rz(2.1329761) q[1];
x q[2];
rz(0.86117427) q[3];
sx q[3];
rz(-2.3434601) q[3];
sx q[3];
rz(2.8482311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5143738) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(-2.880704) q[2];
rz(2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(-2.2055431) q[0];
rz(3.127457) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(0.65151185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053035887) q[0];
sx q[0];
rz(-0.024081973) q[0];
sx q[0];
rz(1.3911029) q[0];
x q[1];
rz(2.3867943) q[2];
sx q[2];
rz(-0.81364606) q[2];
sx q[2];
rz(-0.99572832) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1134909) q[1];
sx q[1];
rz(-0.3416225) q[1];
sx q[1];
rz(-1.2042868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9367847) q[3];
sx q[3];
rz(-1.1567187) q[3];
sx q[3];
rz(0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8462048) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(-2.8335617) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.3778936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4353367) q[0];
sx q[0];
rz(-1.0284817) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3465967) q[2];
sx q[2];
rz(-2.1254052) q[2];
sx q[2];
rz(-2.5650052) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1174406) q[1];
sx q[1];
rz(-1.1785058) q[1];
sx q[1];
rz(-1.8333608) q[1];
rz(0.25791191) q[3];
sx q[3];
rz(-0.57255089) q[3];
sx q[3];
rz(2.8618328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.750981) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(2.509027) q[1];
sx q[1];
rz(-2.0805151) q[1];
sx q[1];
rz(0.17593304) q[1];
rz(2.73545) q[2];
sx q[2];
rz(-2.4293025) q[2];
sx q[2];
rz(0.14355125) q[2];
rz(-0.73344161) q[3];
sx q[3];
rz(-1.0167132) q[3];
sx q[3];
rz(0.50501962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];