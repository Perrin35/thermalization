OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(-0.9932819) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(1.4176523) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7844789) q[0];
sx q[0];
rz(-2.3001101) q[0];
sx q[0];
rz(-0.48805663) q[0];
x q[1];
rz(-2.4383847) q[2];
sx q[2];
rz(-1.1242014) q[2];
sx q[2];
rz(-2.869702) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57342178) q[1];
sx q[1];
rz(-0.69458285) q[1];
sx q[1];
rz(-2.707469) q[1];
x q[2];
rz(-1.3325813) q[3];
sx q[3];
rz(-0.138126) q[3];
sx q[3];
rz(-2.3254243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.950497) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(1.0547868) q[2];
rz(-2.6705006) q[3];
sx q[3];
rz(-1.6867009) q[3];
sx q[3];
rz(-0.86827046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3926113) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(2.9066322) q[0];
rz(1.8670392) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(-0.68769208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0339081) q[0];
sx q[0];
rz(-1.8982197) q[0];
sx q[0];
rz(1.4998687) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1178431) q[2];
sx q[2];
rz(-0.79086429) q[2];
sx q[2];
rz(-0.15463725) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96504922) q[1];
sx q[1];
rz(-1.5582917) q[1];
sx q[1];
rz(3.1360044) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0491203) q[3];
sx q[3];
rz(-2.1824129) q[3];
sx q[3];
rz(3.0870952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4253652) q[2];
sx q[2];
rz(-1.8391515) q[2];
sx q[2];
rz(-0.20935527) q[2];
rz(2.1685205) q[3];
sx q[3];
rz(-2.2416185) q[3];
sx q[3];
rz(-0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3000325) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(2.9969065) q[0];
rz(1.2380098) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(-0.091781052) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84052575) q[0];
sx q[0];
rz(-0.65945461) q[0];
sx q[0];
rz(2.1546138) q[0];
rz(-pi) q[1];
rz(2.1917159) q[2];
sx q[2];
rz(-0.11813049) q[2];
sx q[2];
rz(-1.0850414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8615177) q[1];
sx q[1];
rz(-0.28805486) q[1];
sx q[1];
rz(2.3029598) q[1];
x q[2];
rz(-1.4258521) q[3];
sx q[3];
rz(-2.7053389) q[3];
sx q[3];
rz(0.031883447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7564275) q[2];
sx q[2];
rz(-1.6331208) q[2];
sx q[2];
rz(-2.7623994) q[2];
rz(-0.82342974) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32651153) q[0];
sx q[0];
rz(-2.264475) q[0];
sx q[0];
rz(2.1429578) q[0];
rz(2.6594992) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(1.3179717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40967746) q[0];
sx q[0];
rz(-1.9008667) q[0];
sx q[0];
rz(0.56185742) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58546328) q[2];
sx q[2];
rz(-1.3635067) q[2];
sx q[2];
rz(0.49400615) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6819344) q[1];
sx q[1];
rz(-2.1005957) q[1];
sx q[1];
rz(1.8145241) q[1];
rz(-pi) q[2];
rz(1.0467806) q[3];
sx q[3];
rz(-0.84924997) q[3];
sx q[3];
rz(1.0224354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92129293) q[2];
sx q[2];
rz(-0.90196323) q[2];
sx q[2];
rz(2.656142) q[2];
rz(0.38393936) q[3];
sx q[3];
rz(-1.5860312) q[3];
sx q[3];
rz(0.78401047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-3.0601198) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(-3.1299348) q[0];
rz(-0.13721379) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(1.8395909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0948167) q[0];
sx q[0];
rz(-1.4621549) q[0];
sx q[0];
rz(-2.957445) q[0];
x q[1];
rz(2.1757893) q[2];
sx q[2];
rz(-2.0070672) q[2];
sx q[2];
rz(-2.3406313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1630262) q[1];
sx q[1];
rz(-0.33125505) q[1];
sx q[1];
rz(1.909373) q[1];
rz(-pi) q[2];
rz(0.75826606) q[3];
sx q[3];
rz(-2.6406384) q[3];
sx q[3];
rz(1.4877332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4181218) q[2];
sx q[2];
rz(-1.8401517) q[2];
sx q[2];
rz(0.29329014) q[2];
rz(3.0784741) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(-0.89490923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1191331) q[0];
sx q[0];
rz(-0.28894153) q[0];
sx q[0];
rz(3.1030848) q[0];
rz(1.0844082) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(-0.96673059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16822505) q[0];
sx q[0];
rz(-1.7666139) q[0];
sx q[0];
rz(-1.7652579) q[0];
rz(-3.128563) q[2];
sx q[2];
rz(-1.8378864) q[2];
sx q[2];
rz(-2.992127) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37361426) q[1];
sx q[1];
rz(-1.9100338) q[1];
sx q[1];
rz(-1.148073) q[1];
rz(-1.8943664) q[3];
sx q[3];
rz(-2.4651291) q[3];
sx q[3];
rz(-2.885779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4806369) q[2];
sx q[2];
rz(-0.2350685) q[2];
sx q[2];
rz(-2.1601775) q[2];
rz(1.4977945) q[3];
sx q[3];
rz(-1.2621597) q[3];
sx q[3];
rz(-1.5952544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82866955) q[0];
sx q[0];
rz(-0.68083119) q[0];
sx q[0];
rz(-2.8269826) q[0];
rz(0.18868748) q[1];
sx q[1];
rz(-2.7068832) q[1];
sx q[1];
rz(2.6729118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9148736) q[0];
sx q[0];
rz(-1.3283936) q[0];
sx q[0];
rz(-1.4366494) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5596703) q[2];
sx q[2];
rz(-2.9738148) q[2];
sx q[2];
rz(1.1388701) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78476731) q[1];
sx q[1];
rz(-1.2677437) q[1];
sx q[1];
rz(2.3222938) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46814274) q[3];
sx q[3];
rz(-1.2791514) q[3];
sx q[3];
rz(2.1848752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20588747) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(2.1708798) q[2];
rz(-0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(0.71603388) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69990528) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(1.4458789) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(-3.0605002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4202284) q[0];
sx q[0];
rz(-2.2330771) q[0];
sx q[0];
rz(-1.2032979) q[0];
rz(2.4862635) q[2];
sx q[2];
rz(-1.809279) q[2];
sx q[2];
rz(0.9376038) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9633741) q[1];
sx q[1];
rz(-1.4556754) q[1];
sx q[1];
rz(-1.2769775) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.567566) q[3];
sx q[3];
rz(-1.9732631) q[3];
sx q[3];
rz(-0.57202475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1991835) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(-0.24277631) q[2];
rz(-1.5822004) q[3];
sx q[3];
rz(-2.715761) q[3];
sx q[3];
rz(1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9026069) q[0];
sx q[0];
rz(-1.0317529) q[0];
sx q[0];
rz(1.8865939) q[0];
rz(1.7651419) q[1];
sx q[1];
rz(-1.0245198) q[1];
sx q[1];
rz(-0.66744101) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4723929) q[0];
sx q[0];
rz(-2.071864) q[0];
sx q[0];
rz(0.38031148) q[0];
rz(1.5282643) q[2];
sx q[2];
rz(-2.0669524) q[2];
sx q[2];
rz(0.68788487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63782843) q[1];
sx q[1];
rz(-2.1385647) q[1];
sx q[1];
rz(2.7953887) q[1];
x q[2];
rz(1.8967751) q[3];
sx q[3];
rz(-1.2080844) q[3];
sx q[3];
rz(-2.6397702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5218375) q[2];
sx q[2];
rz(-0.73885584) q[2];
sx q[2];
rz(2.6832306) q[2];
rz(1.6058263) q[3];
sx q[3];
rz(-1.0777487) q[3];
sx q[3];
rz(2.2569236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430724) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(2.7761053) q[0];
rz(-0.99994031) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(1.4568636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92577584) q[0];
sx q[0];
rz(-2.121637) q[0];
sx q[0];
rz(2.9047545) q[0];
rz(-pi) q[1];
rz(-1.959895) q[2];
sx q[2];
rz(-0.28841296) q[2];
sx q[2];
rz(2.4450977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0914336) q[1];
sx q[1];
rz(-0.83178565) q[1];
sx q[1];
rz(1.8065401) q[1];
rz(-pi) q[2];
rz(2.6014464) q[3];
sx q[3];
rz(-2.3688101) q[3];
sx q[3];
rz(2.5908822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.140787) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(-2.136039) q[2];
rz(2.4908861) q[3];
sx q[3];
rz(-1.7020099) q[3];
sx q[3];
rz(-2.2910291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093813048) q[0];
sx q[0];
rz(-0.78421264) q[0];
sx q[0];
rz(-1.0017851) q[0];
rz(-1.4108989) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(1.9000353) q[2];
sx q[2];
rz(-1.4598979) q[2];
sx q[2];
rz(2.0171063) q[2];
rz(1.4776354) q[3];
sx q[3];
rz(-1.6747083) q[3];
sx q[3];
rz(-1.2887521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
