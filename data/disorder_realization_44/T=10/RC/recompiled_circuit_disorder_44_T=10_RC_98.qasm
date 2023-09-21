OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44088988) q[0];
sx q[0];
rz(-2.7452677) q[0];
sx q[0];
rz(-1.8519782) q[0];
rz(-2.2719703) q[2];
sx q[2];
rz(-0.50719559) q[2];
sx q[2];
rz(1.6091572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9143608) q[1];
sx q[1];
rz(-2.025361) q[1];
sx q[1];
rz(0.55959065) q[1];
rz(2.3832541) q[3];
sx q[3];
rz(-1.3876649) q[3];
sx q[3];
rz(1.5116215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-3.120378) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(0.83591998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056136925) q[0];
sx q[0];
rz(-1.5526999) q[0];
sx q[0];
rz(1.4319112) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40936562) q[2];
sx q[2];
rz(-1.0592807) q[2];
sx q[2];
rz(0.53586938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9906209) q[1];
sx q[1];
rz(-0.70578209) q[1];
sx q[1];
rz(-1.1433931) q[1];
x q[2];
rz(2.7153035) q[3];
sx q[3];
rz(-2.3549035) q[3];
sx q[3];
rz(-2.8601437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(0.35955444) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-0.4371117) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1431883) q[0];
sx q[0];
rz(-0.30448738) q[0];
sx q[0];
rz(1.8633153) q[0];
rz(1.1630467) q[2];
sx q[2];
rz(-1.7951269) q[2];
sx q[2];
rz(0.11220223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5324085) q[1];
sx q[1];
rz(-0.44645616) q[1];
sx q[1];
rz(2.8162454) q[1];
rz(-pi) q[2];
rz(1.0760355) q[3];
sx q[3];
rz(-1.5177739) q[3];
sx q[3];
rz(0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(-1.0220698) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(2.7220272) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(-3.006382) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-2.9503126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0606196) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(0.087555126) q[0];
x q[1];
rz(-0.27387597) q[2];
sx q[2];
rz(-2.2857776) q[2];
sx q[2];
rz(0.82211923) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7820218) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(2.3624079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6077765) q[3];
sx q[3];
rz(-1.8752408) q[3];
sx q[3];
rz(0.75418762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-1.0106687) q[2];
rz(-0.7615532) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(2.9060569) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199729) q[0];
sx q[0];
rz(-3.0191506) q[0];
sx q[0];
rz(-0.58453154) q[0];
rz(-pi) q[1];
rz(0.82053484) q[2];
sx q[2];
rz(-1.8168601) q[2];
sx q[2];
rz(0.29124242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3481808) q[1];
sx q[1];
rz(-2.9017397) q[1];
sx q[1];
rz(-2.1972149) q[1];
x q[2];
rz(2.186741) q[3];
sx q[3];
rz(-1.5408437) q[3];
sx q[3];
rz(-1.1101462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(1.3657773) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(-2.1877066) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(2.0772207) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-0.37429601) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74852809) q[0];
sx q[0];
rz(-0.52951282) q[0];
sx q[0];
rz(2.4783496) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49403814) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(-1.5766694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6090138) q[1];
sx q[1];
rz(-0.62090579) q[1];
sx q[1];
rz(-0.46355526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80769844) q[3];
sx q[3];
rz(-0.85113111) q[3];
sx q[3];
rz(0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8950243) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(0.95823112) q[2];
rz(2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(-1.3100756) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3416672) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-3.0068586) q[0];
rz(-pi) q[1];
rz(1.9016978) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(-0.14724018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.29855) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(0.33903867) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55540107) q[3];
sx q[3];
rz(-1.5474833) q[3];
sx q[3];
rz(-2.5454552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16335547) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(1.3423052) q[0];
rz(2.2539027) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(-2.2788252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0975768) q[1];
sx q[1];
rz(-2.3313287) q[1];
sx q[1];
rz(0.047659831) q[1];
rz(-pi) q[2];
rz(-1.3109342) q[3];
sx q[3];
rz(-1.6901008) q[3];
sx q[3];
rz(2.3343149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(0.35593629) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5091771) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-0.40503043) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.2776432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67280806) q[0];
sx q[0];
rz(-2.0622258) q[0];
sx q[0];
rz(1.3093033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0870886) q[2];
sx q[2];
rz(-0.93855575) q[2];
sx q[2];
rz(1.5916923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35754044) q[1];
sx q[1];
rz(-1.5338384) q[1];
sx q[1];
rz(1.1251015) q[1];
x q[2];
rz(1.1489264) q[3];
sx q[3];
rz(-0.78242362) q[3];
sx q[3];
rz(-2.9187834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3383125) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(-1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(2.0675802) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-3.0293368) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(2.9311438) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1505517) q[0];
sx q[0];
rz(-1.7518839) q[0];
sx q[0];
rz(-2.6291356) q[0];
rz(2.9774882) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(1.1526398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75281843) q[1];
sx q[1];
rz(-0.37467271) q[1];
sx q[1];
rz(0.59605662) q[1];
x q[2];
rz(-0.75746234) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(-2.7528742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(2.9343228) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5719941) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(0.81644425) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(-1.2522092) q[2];
sx q[2];
rz(-0.49124419) q[2];
sx q[2];
rz(-1.0791525) q[2];
rz(1.1917226) q[3];
sx q[3];
rz(-0.76615292) q[3];
sx q[3];
rz(-2.5828008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
