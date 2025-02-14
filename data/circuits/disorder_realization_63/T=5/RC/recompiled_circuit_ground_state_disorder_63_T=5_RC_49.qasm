OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(1.0898606) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(-0.90322948) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2555465) q[0];
sx q[0];
rz(-2.1849792) q[0];
sx q[0];
rz(-0.79844676) q[0];
rz(1.5527234) q[2];
sx q[2];
rz(-0.68328349) q[2];
sx q[2];
rz(-2.9073213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.049202327) q[1];
sx q[1];
rz(-2.4396439) q[1];
sx q[1];
rz(-2.5939221) q[1];
rz(-1.0711477) q[3];
sx q[3];
rz(-1.6766492) q[3];
sx q[3];
rz(-1.203804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1746154) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(1.2906661) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(2.2490833) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21890579) q[0];
sx q[0];
rz(-1.7480047) q[0];
sx q[0];
rz(-0.025434504) q[0];
x q[1];
rz(-2.4878534) q[2];
sx q[2];
rz(-1.8934665) q[2];
sx q[2];
rz(1.309909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.635517) q[1];
sx q[1];
rz(-0.92787023) q[1];
sx q[1];
rz(-1.9159813) q[1];
rz(1.1970911) q[3];
sx q[3];
rz(-1.804816) q[3];
sx q[3];
rz(-1.6683874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9925925) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(-2.0089669) q[2];
rz(-2.4752786) q[3];
sx q[3];
rz(-1.7814813) q[3];
sx q[3];
rz(-1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76729524) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(0.98180109) q[0];
rz(2.1994195) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(2.2312677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3142598) q[0];
sx q[0];
rz(-1.1480486) q[0];
sx q[0];
rz(0.5323005) q[0];
x q[1];
rz(-1.4508171) q[2];
sx q[2];
rz(-1.8455122) q[2];
sx q[2];
rz(-3.0697696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83214789) q[1];
sx q[1];
rz(-0.4570241) q[1];
sx q[1];
rz(0.8330658) q[1];
rz(-pi) q[2];
rz(-1.1317838) q[3];
sx q[3];
rz(-1.7903504) q[3];
sx q[3];
rz(1.1846402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2008449) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(-2.9546837) q[2];
rz(0.16380353) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(2.4439268) q[0];
rz(0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.9465416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89681399) q[0];
sx q[0];
rz(-1.91511) q[0];
sx q[0];
rz(-2.1953775) q[0];
rz(-2.5279453) q[2];
sx q[2];
rz(-0.23229182) q[2];
sx q[2];
rz(-2.7181546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4711485) q[1];
sx q[1];
rz(-2.9440144) q[1];
sx q[1];
rz(-1.0001282) q[1];
rz(3.0680858) q[3];
sx q[3];
rz(-1.2871082) q[3];
sx q[3];
rz(0.43108854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67808548) q[2];
sx q[2];
rz(-1.4231851) q[2];
sx q[2];
rz(-0.19742337) q[2];
rz(-0.10489634) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(0.088002861) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(1.0092258) q[0];
rz(3.1127473) q[1];
sx q[1];
rz(-1.6639158) q[1];
sx q[1];
rz(-2.4829594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742666) q[0];
sx q[0];
rz(-1.9443041) q[0];
sx q[0];
rz(-2.0167355) q[0];
rz(-0.20435996) q[2];
sx q[2];
rz(-2.0484701) q[2];
sx q[2];
rz(2.7088239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11396961) q[1];
sx q[1];
rz(-1.4024892) q[1];
sx q[1];
rz(-1.6644003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3405209) q[3];
sx q[3];
rz(-2.252823) q[3];
sx q[3];
rz(-1.7868228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.37904) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(-2.9075882) q[2];
rz(2.0512569) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91319084) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(-1.7826805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8131065) q[0];
sx q[0];
rz(-0.66603249) q[0];
sx q[0];
rz(3.0124927) q[0];
x q[1];
rz(-1.3471021) q[2];
sx q[2];
rz(-1.3094721) q[2];
sx q[2];
rz(1.4841207) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7715367) q[1];
sx q[1];
rz(-3.1395457) q[1];
sx q[1];
rz(-2.5084346) q[1];
rz(-pi) q[2];
rz(0.5868191) q[3];
sx q[3];
rz(-0.53884655) q[3];
sx q[3];
rz(0.19308819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4141554) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(-2.9690913) q[2];
rz(-2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(0.33997047) q[0];
rz(-2.8151457) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(1.9065769) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3351072) q[0];
sx q[0];
rz(-2.3705784) q[0];
sx q[0];
rz(-0.54698995) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9307617) q[2];
sx q[2];
rz(-2.0034997) q[2];
sx q[2];
rz(2.3683647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82385027) q[1];
sx q[1];
rz(-0.47951937) q[1];
sx q[1];
rz(3.13003) q[1];
x q[2];
rz(0.39374473) q[3];
sx q[3];
rz(-1.7753768) q[3];
sx q[3];
rz(-1.0130628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1116703) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(2.663747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70917201) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(1.1093371) q[0];
rz(2.0093911) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(2.1999377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3170211) q[0];
sx q[0];
rz(-2.5618658) q[0];
sx q[0];
rz(-1.8648022) q[0];
rz(-pi) q[1];
rz(1.2852816) q[2];
sx q[2];
rz(-1.4219577) q[2];
sx q[2];
rz(-0.51644737) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1600765) q[1];
sx q[1];
rz(-1.4443411) q[1];
sx q[1];
rz(-0.92595788) q[1];
rz(1.3579943) q[3];
sx q[3];
rz(-1.1929242) q[3];
sx q[3];
rz(1.6869643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7534916) q[2];
sx q[2];
rz(-2.1145623) q[2];
sx q[2];
rz(-1.5388185) q[2];
rz(1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(-2.1737461) q[0];
rz(-1.2713426) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(3.079792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4000452) q[0];
sx q[0];
rz(-1.9572659) q[0];
sx q[0];
rz(2.5611112) q[0];
rz(-0.25904663) q[2];
sx q[2];
rz(-1.3767837) q[2];
sx q[2];
rz(1.2736748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3810515) q[1];
sx q[1];
rz(-0.88812056) q[1];
sx q[1];
rz(1.6939916) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7155719) q[3];
sx q[3];
rz(-0.4288097) q[3];
sx q[3];
rz(-2.9675067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1900078) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(-3.0657213) q[2];
rz(1.1672945) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(-2.8372138) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17700125) q[0];
sx q[0];
rz(-0.27934203) q[0];
sx q[0];
rz(-2.8569073) q[0];
rz(0.074706569) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(-1.7237192) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3475932) q[0];
sx q[0];
rz(-2.3132304) q[0];
sx q[0];
rz(1.30778) q[0];
rz(-2.9921164) q[2];
sx q[2];
rz(-1.0188531) q[2];
sx q[2];
rz(1.1750482) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7602882) q[1];
sx q[1];
rz(-2.3592296) q[1];
sx q[1];
rz(-2.5336877) q[1];
rz(-pi) q[2];
rz(0.033662658) q[3];
sx q[3];
rz(-0.92516781) q[3];
sx q[3];
rz(1.0317623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8031926) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(-1.2791862) q[2];
rz(-2.159436) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(-2.2568683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2411156) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(2.1263532) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(1.8335291) q[2];
sx q[2];
rz(-1.4663525) q[2];
sx q[2];
rz(2.732085) q[2];
rz(-2.1463263) q[3];
sx q[3];
rz(-2.223816) q[3];
sx q[3];
rz(1.73903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
