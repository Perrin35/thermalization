OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0187674) q[0];
sx q[0];
rz(-3.1376165) q[0];
sx q[0];
rz(3.0206326) q[0];
rz(0.51789415) q[1];
sx q[1];
rz(-0.33101141) q[1];
sx q[1];
rz(1.9537227) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56138229) q[0];
sx q[0];
rz(-2.8299709) q[0];
sx q[0];
rz(2.8557957) q[0];
rz(1.0446494) q[2];
sx q[2];
rz(-1.6974715) q[2];
sx q[2];
rz(-2.130545) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49379865) q[1];
sx q[1];
rz(-1.0435074) q[1];
sx q[1];
rz(-2.6721558) q[1];
rz(-pi) q[2];
rz(2.5099947) q[3];
sx q[3];
rz(-1.609841) q[3];
sx q[3];
rz(-1.6717048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0594222) q[2];
sx q[2];
rz(-2.1990081) q[2];
sx q[2];
rz(1.7131294) q[2];
rz(1.3905455) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(-2.9318504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21935894) q[0];
sx q[0];
rz(-1.5606422) q[0];
sx q[0];
rz(-1.4571762) q[0];
rz(-2.4600785) q[1];
sx q[1];
rz(-0.69029713) q[1];
sx q[1];
rz(-0.95726475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663366) q[0];
sx q[0];
rz(-2.0342376) q[0];
sx q[0];
rz(-1.8011679) q[0];
x q[1];
rz(0.817746) q[2];
sx q[2];
rz(-1.3829573) q[2];
sx q[2];
rz(-0.85603324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71979501) q[1];
sx q[1];
rz(-1.52175) q[1];
sx q[1];
rz(1.979046) q[1];
rz(-pi) q[2];
rz(0.010250895) q[3];
sx q[3];
rz(-1.6396549) q[3];
sx q[3];
rz(2.0673925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7822781) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(0.90685833) q[2];
rz(-0.67306972) q[3];
sx q[3];
rz(-2.7570717) q[3];
sx q[3];
rz(-1.2955906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(0.54687706) q[0];
rz(-1.3705672) q[1];
sx q[1];
rz(-1.9799045) q[1];
sx q[1];
rz(-3.0414157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42368314) q[0];
sx q[0];
rz(-0.29742213) q[0];
sx q[0];
rz(-1.1017767) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9630394) q[2];
sx q[2];
rz(-1.0970976) q[2];
sx q[2];
rz(0.77609962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5838769) q[1];
sx q[1];
rz(-2.2889232) q[1];
sx q[1];
rz(1.2664944) q[1];
x q[2];
rz(0.46296316) q[3];
sx q[3];
rz(-1.1826666) q[3];
sx q[3];
rz(-2.8007647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3935516) q[2];
sx q[2];
rz(-1.9358205) q[2];
sx q[2];
rz(-2.098341) q[2];
rz(-0.68759632) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461175) q[0];
sx q[0];
rz(-2.6961374) q[0];
sx q[0];
rz(-1.0365781) q[0];
rz(-1.2713894) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(-1.8545256) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66726738) q[0];
sx q[0];
rz(-2.0653353) q[0];
sx q[0];
rz(-0.70566341) q[0];
rz(2.4802568) q[2];
sx q[2];
rz(-1.9171538) q[2];
sx q[2];
rz(2.6544184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8657664) q[1];
sx q[1];
rz(-0.681502) q[1];
sx q[1];
rz(1.1033415) q[1];
x q[2];
rz(-0.98388328) q[3];
sx q[3];
rz(-0.5183087) q[3];
sx q[3];
rz(0.89215088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.53671545) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(-0.019901179) q[2];
rz(-2.5455635) q[3];
sx q[3];
rz(-2.8607131) q[3];
sx q[3];
rz(-1.178044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48536456) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(1.9670991) q[0];
rz(-0.30612048) q[1];
sx q[1];
rz(-2.0968292) q[1];
sx q[1];
rz(-1.3694481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58254439) q[0];
sx q[0];
rz(-1.567548) q[0];
sx q[0];
rz(1.5772468) q[0];
x q[1];
rz(-2.1691983) q[2];
sx q[2];
rz(-1.4273424) q[2];
sx q[2];
rz(0.17696887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7983129) q[1];
sx q[1];
rz(-1.7213744) q[1];
sx q[1];
rz(1.2231138) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6275614) q[3];
sx q[3];
rz(-1.2054878) q[3];
sx q[3];
rz(-1.1082054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29641637) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(0.80294341) q[2];
rz(-2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1882316) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(-0.11235919) q[0];
rz(0.27680963) q[1];
sx q[1];
rz(-1.9748297) q[1];
sx q[1];
rz(2.4940431) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1153664) q[0];
sx q[0];
rz(-2.1121968) q[0];
sx q[0];
rz(0.39570915) q[0];
rz(-pi) q[1];
rz(2.4729257) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(2.6636555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4989483) q[1];
sx q[1];
rz(-2.0181839) q[1];
sx q[1];
rz(-1.7018464) q[1];
x q[2];
rz(1.8065436) q[3];
sx q[3];
rz(-1.4000579) q[3];
sx q[3];
rz(2.6125107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67927805) q[2];
sx q[2];
rz(-0.15028149) q[2];
sx q[2];
rz(1.798299) q[2];
rz(2.6278833) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(-0.59521365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419471) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(0.42022002) q[0];
rz(1.7814024) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(-2.3374048) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4046832) q[0];
sx q[0];
rz(-2.1054287) q[0];
sx q[0];
rz(-0.44628365) q[0];
rz(-pi) q[1];
rz(-0.28557318) q[2];
sx q[2];
rz(-0.80508666) q[2];
sx q[2];
rz(-2.7325163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7090164) q[1];
sx q[1];
rz(-0.99919277) q[1];
sx q[1];
rz(2.9663621) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4057833) q[3];
sx q[3];
rz(-1.4632153) q[3];
sx q[3];
rz(-3.0522814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5516025) q[2];
sx q[2];
rz(-2.1253822) q[2];
sx q[2];
rz(-1.972398) q[2];
rz(0.17360887) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(-1.199523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6731897) q[0];
sx q[0];
rz(-2.115374) q[0];
sx q[0];
rz(-1.9386559) q[0];
rz(2.9007593) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(-0.9333207) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627181) q[0];
sx q[0];
rz(-1.2347455) q[0];
sx q[0];
rz(-1.1707992) q[0];
x q[1];
rz(-1.9674932) q[2];
sx q[2];
rz(-2.5317414) q[2];
sx q[2];
rz(1.6624474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9874064) q[1];
sx q[1];
rz(-1.1917104) q[1];
sx q[1];
rz(1.3554708) q[1];
rz(-pi) q[2];
rz(-0.15214129) q[3];
sx q[3];
rz(-1.0120494) q[3];
sx q[3];
rz(1.5882176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11747083) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(-0.13548279) q[2];
rz(-0.99902117) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(-2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713772) q[0];
sx q[0];
rz(-2.0129634) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(2.0595835) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(1.5622004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45856548) q[0];
sx q[0];
rz(-0.54046721) q[0];
sx q[0];
rz(-1.6349273) q[0];
x q[1];
rz(0.57317971) q[2];
sx q[2];
rz(-0.64116353) q[2];
sx q[2];
rz(1.6784843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.434103) q[1];
sx q[1];
rz(-1.130778) q[1];
sx q[1];
rz(-0.36807366) q[1];
rz(2.2213577) q[3];
sx q[3];
rz(-1.2830858) q[3];
sx q[3];
rz(-0.77863151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9595327) q[2];
sx q[2];
rz(-1.5675194) q[2];
sx q[2];
rz(2.5173397) q[2];
rz(-0.67638451) q[3];
sx q[3];
rz(-1.9703777) q[3];
sx q[3];
rz(-0.83522767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4312209) q[0];
sx q[0];
rz(-1.8932764) q[0];
sx q[0];
rz(-1.7927908) q[0];
rz(2.8377332) q[1];
sx q[1];
rz(-1.5307348) q[1];
sx q[1];
rz(2.2644728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6346587) q[0];
sx q[0];
rz(-1.7153694) q[0];
sx q[0];
rz(0.73214447) q[0];
rz(1.2563324) q[2];
sx q[2];
rz(-2.8012677) q[2];
sx q[2];
rz(0.2257502) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7925229) q[1];
sx q[1];
rz(-1.190509) q[1];
sx q[1];
rz(-2.040634) q[1];
rz(0.83721807) q[3];
sx q[3];
rz(-2.3851237) q[3];
sx q[3];
rz(-2.9314007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0554241) q[2];
sx q[2];
rz(-0.40579) q[2];
sx q[2];
rz(-2.1093192) q[2];
rz(-1.0885193) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(2.3672339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9536204) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(-1.8232518) q[1];
sx q[1];
rz(-1.4069469) q[1];
sx q[1];
rz(-2.5896752) q[1];
rz(-1.8404519) q[2];
sx q[2];
rz(-1.1674623) q[2];
sx q[2];
rz(1.8415378) q[2];
rz(-1.134675) q[3];
sx q[3];
rz(-1.5584617) q[3];
sx q[3];
rz(-0.8212318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
