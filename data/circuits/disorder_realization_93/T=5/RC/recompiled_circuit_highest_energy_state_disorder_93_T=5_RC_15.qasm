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
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(-1.4135345) q[0];
rz(0.36355525) q[1];
sx q[1];
rz(-1.3132562) q[1];
sx q[1];
rz(-0.29120905) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.04714) q[0];
sx q[0];
rz(-1.5507076) q[0];
sx q[0];
rz(-1.8909341) q[0];
rz(1.232336) q[2];
sx q[2];
rz(-0.3582274) q[2];
sx q[2];
rz(-3.0436277) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8785581) q[1];
sx q[1];
rz(-0.59459762) q[1];
sx q[1];
rz(-2.8697398) q[1];
rz(-2.324027) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(2.7369902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0061965813) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(-0.023999365) q[2];
rz(0.35605797) q[3];
sx q[3];
rz(-0.67432109) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085670797) q[0];
sx q[0];
rz(-0.69127685) q[0];
sx q[0];
rz(2.5066277) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(-2.2244804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0046154) q[0];
sx q[0];
rz(-0.03923035) q[0];
sx q[0];
rz(-1.0835539) q[0];
x q[1];
rz(-1.00441) q[2];
sx q[2];
rz(-0.76493636) q[2];
sx q[2];
rz(1.9435918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0113127) q[1];
sx q[1];
rz(-0.44751274) q[1];
sx q[1];
rz(-1.7496197) q[1];
rz(-0.95083745) q[3];
sx q[3];
rz(-1.3966433) q[3];
sx q[3];
rz(0.0096498409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94747535) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(-2.4083162) q[2];
rz(-2.0745847) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5697524) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(2.6285062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.899174) q[0];
sx q[0];
rz(-2.0089825) q[0];
sx q[0];
rz(0.04695453) q[0];
rz(1.1703067) q[2];
sx q[2];
rz(-2.4520564) q[2];
sx q[2];
rz(-2.6282118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9218623) q[1];
sx q[1];
rz(-0.37223909) q[1];
sx q[1];
rz(-1.2138741) q[1];
rz(-pi) q[2];
rz(0.49969466) q[3];
sx q[3];
rz(-1.5958188) q[3];
sx q[3];
rz(2.6270888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1055866) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(-0.038724381) q[2];
rz(-2.683908) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(-1.800644) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(-2.0204954) q[0];
rz(-0.7577678) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(-2.4574492) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0476889) q[0];
sx q[0];
rz(-2.2950116) q[0];
sx q[0];
rz(-2.8709477) q[0];
rz(-0.010104601) q[2];
sx q[2];
rz(-1.5636383) q[2];
sx q[2];
rz(0.3601851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18312081) q[1];
sx q[1];
rz(-2.712912) q[1];
sx q[1];
rz(-0.76629345) q[1];
rz(-2.6287967) q[3];
sx q[3];
rz(-1.0362451) q[3];
sx q[3];
rz(-3.0844546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(2.8810775) q[2];
rz(0.83812964) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-0.014613541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218006) q[0];
sx q[0];
rz(-1.7638693) q[0];
sx q[0];
rz(-2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.1621954) q[1];
sx q[1];
rz(-1.6872663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5430253) q[0];
sx q[0];
rz(-1.784458) q[0];
sx q[0];
rz(1.5419307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32938214) q[2];
sx q[2];
rz(-2.9940146) q[2];
sx q[2];
rz(-1.3775228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37794681) q[1];
sx q[1];
rz(-1.0275405) q[1];
sx q[1];
rz(1.8733979) q[1];
rz(-pi) q[2];
rz(-0.14438676) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(-0.99100366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55383596) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(-2.7510263) q[2];
rz(-2.8824814) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(-0.34671569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(-0.64467347) q[0];
rz(1.3020172) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(2.7929746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020418305) q[0];
sx q[0];
rz(-2.0266337) q[0];
sx q[0];
rz(1.0555847) q[0];
x q[1];
rz(0.63747823) q[2];
sx q[2];
rz(-2.5036252) q[2];
sx q[2];
rz(-1.5359985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8585414) q[1];
sx q[1];
rz(-1.9894454) q[1];
sx q[1];
rz(-2.8178697) q[1];
rz(-pi) q[2];
rz(-2.8556267) q[3];
sx q[3];
rz(-0.12187258) q[3];
sx q[3];
rz(-1.414145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0384486) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(0.038334282) q[2];
rz(-2.1740055) q[3];
sx q[3];
rz(-1.4097694) q[3];
sx q[3];
rz(-2.7819395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(-1.7346802) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(0.73307347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382815) q[0];
sx q[0];
rz(-1.065863) q[0];
sx q[0];
rz(2.5236593) q[0];
x q[1];
rz(2.2233811) q[2];
sx q[2];
rz(-0.84647734) q[2];
sx q[2];
rz(-1.0670964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.275707) q[1];
sx q[1];
rz(-1.2162677) q[1];
sx q[1];
rz(-1.2178631) q[1];
rz(-pi) q[2];
x q[2];
rz(1.319995) q[3];
sx q[3];
rz(-1.7989252) q[3];
sx q[3];
rz(-3.1056014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1072032) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(-0.011073152) q[2];
rz(0.30558807) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(-1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333862) q[0];
sx q[0];
rz(-2.5337063) q[0];
sx q[0];
rz(2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(0.098800585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042511333) q[0];
sx q[0];
rz(-2.3020199) q[0];
sx q[0];
rz(2.4057516) q[0];
rz(-pi) q[1];
rz(-2.4315028) q[2];
sx q[2];
rz(-2.3259893) q[2];
sx q[2];
rz(2.2888773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0055726) q[1];
sx q[1];
rz(-0.51678777) q[1];
sx q[1];
rz(-2.6173475) q[1];
rz(-3.0116073) q[3];
sx q[3];
rz(-1.49125) q[3];
sx q[3];
rz(-1.4628439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8431479) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(2.0419545) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(-2.5449424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1140887) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(-1.9644568) q[0];
rz(-1.6674532) q[1];
sx q[1];
rz(-2.2087704) q[1];
sx q[1];
rz(1.6966381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5367573) q[0];
sx q[0];
rz(-2.6492181) q[0];
sx q[0];
rz(-0.62712609) q[0];
rz(-pi) q[1];
rz(-0.58622245) q[2];
sx q[2];
rz(-2.2811523) q[2];
sx q[2];
rz(1.1733871) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7921389) q[1];
sx q[1];
rz(-1.0568403) q[1];
sx q[1];
rz(-3.0626014) q[1];
rz(-pi) q[2];
rz(-1.2304495) q[3];
sx q[3];
rz(-2.4370464) q[3];
sx q[3];
rz(-1.6753538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(0.17189279) q[2];
rz(-2.3520172) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(1.3353698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563141) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(0.54157448) q[0];
rz(1.9821292) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(2.3084739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492958) q[0];
sx q[0];
rz(-0.98092666) q[0];
sx q[0];
rz(-1.2372531) q[0];
rz(-pi) q[1];
rz(1.2255185) q[2];
sx q[2];
rz(-1.0547027) q[2];
sx q[2];
rz(2.6654625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45634746) q[1];
sx q[1];
rz(-2.4083071) q[1];
sx q[1];
rz(-1.228778) q[1];
rz(-pi) q[2];
rz(2.2984357) q[3];
sx q[3];
rz(-1.572248) q[3];
sx q[3];
rz(2.7160318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.021576015) q[2];
sx q[2];
rz(-2.3638201) q[2];
sx q[2];
rz(-2.1982101) q[2];
rz(2.0639482) q[3];
sx q[3];
rz(-1.5871983) q[3];
sx q[3];
rz(2.8130892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42644603) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(0.23964755) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(2.2108261) q[2];
sx q[2];
rz(-0.98251828) q[2];
sx q[2];
rz(0.3457014) q[2];
rz(-0.31044337) q[3];
sx q[3];
rz(-1.3287928) q[3];
sx q[3];
rz(1.2850858) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
