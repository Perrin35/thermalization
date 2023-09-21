OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.278468) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(-1.6186884) q[0];
rz(-2.0179022) q[2];
sx q[2];
rz(-1.2667709) q[2];
sx q[2];
rz(2.6927039) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(-1.0784472) q[1];
x q[2];
rz(-0.11043926) q[3];
sx q[3];
rz(-0.38023284) q[3];
sx q[3];
rz(-1.5649753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.3551691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38793135) q[0];
sx q[0];
rz(-0.61457115) q[0];
sx q[0];
rz(1.4580926) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16263527) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.6441117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5835186) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(3.0706057) q[1];
rz(-pi) q[2];
rz(-1.1739028) q[3];
sx q[3];
rz(-2.9429884) q[3];
sx q[3];
rz(1.393115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(-2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(1.0650939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168034) q[0];
sx q[0];
rz(-1.5771126) q[0];
sx q[0];
rz(1.2114695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2421354) q[2];
sx q[2];
rz(-2.3231635) q[2];
sx q[2];
rz(-0.98086548) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0703735) q[1];
sx q[1];
rz(-2.1871236) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5866833) q[0];
sx q[0];
rz(-2.0218098) q[0];
sx q[0];
rz(0.45278544) q[0];
x q[1];
rz(-0.63979062) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(1.7196136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0032469) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(-2.7702615) q[1];
rz(-pi) q[2];
rz(-2.6075881) q[3];
sx q[3];
rz(-2.7760193) q[3];
sx q[3];
rz(3.0766069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0304612) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.8146851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6700232) q[0];
sx q[0];
rz(-1.4172535) q[0];
sx q[0];
rz(-2.3877386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97557108) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(3.1095568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8618968) q[1];
sx q[1];
rz(-1.3394636) q[1];
sx q[1];
rz(-2.8054603) q[1];
x q[2];
rz(-2.1760686) q[3];
sx q[3];
rz(-2*pi/13) q[3];
sx q[3];
rz(2.0616639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-2.025827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9003446) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(-2.3099398) q[0];
rz(-pi) q[1];
rz(-2.240928) q[2];
sx q[2];
rz(-1.2256983) q[2];
sx q[2];
rz(0.64054856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(0.58847217) q[1];
rz(-1.3917771) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(0.37187809) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-2.5792714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092098504) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(-0.90256079) q[0];
x q[1];
rz(-1.6717031) q[2];
sx q[2];
rz(-1.2404053) q[2];
sx q[2];
rz(2.908387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.454969) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(1.4229694) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.087248487) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.5213535) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884468) q[0];
sx q[0];
rz(-1.7674812) q[0];
sx q[0];
rz(-1.6403857) q[0];
rz(-pi) q[1];
rz(-1.0993768) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(2.9277756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(0.08463879) q[1];
x q[2];
rz(3.0275214) q[3];
sx q[3];
rz(-2.2469006) q[3];
sx q[3];
rz(-1.8166208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(1.1307905) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(0.27063453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(2.4256698) q[0];
x q[1];
rz(-0.46220772) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(2.8414937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2808025) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.526282) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4326914) q[3];
sx q[3];
rz(-2.0376251) q[3];
sx q[3];
rz(0.69463581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4821645) q[0];
sx q[0];
rz(-2.6539301) q[0];
sx q[0];
rz(-1.6517261) q[0];
x q[1];
rz(-1.5779183) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(0.42051007) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(0.93025031) q[1];
rz(-1.6845735) q[3];
sx q[3];
rz(-0.69996951) q[3];
sx q[3];
rz(1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(-2.8137394) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-0.83256759) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(1.238263) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(2.3776688) q[3];
sx q[3];
rz(-0.36119701) q[3];
sx q[3];
rz(1.2154538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
