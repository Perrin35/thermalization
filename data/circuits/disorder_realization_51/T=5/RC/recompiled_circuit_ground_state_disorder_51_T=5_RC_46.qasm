OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12282523) q[0];
sx q[0];
rz(3.1376165) q[0];
sx q[0];
rz(9.3038179) q[0];
rz(-2.6236985) q[1];
sx q[1];
rz(-2.8105812) q[1];
sx q[1];
rz(-1.9537227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26195458) q[0];
sx q[0];
rz(-1.2722135) q[0];
sx q[0];
rz(-1.4802329) q[0];
rz(-2.9953674) q[2];
sx q[2];
rz(-2.0922961) q[2];
sx q[2];
rz(2.5086049) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.647794) q[1];
sx q[1];
rz(-2.0980853) q[1];
sx q[1];
rz(0.46943689) q[1];
rz(1.5224321) q[3];
sx q[3];
rz(-0.93975583) q[3];
sx q[3];
rz(-0.12945718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0594222) q[2];
sx q[2];
rz(-2.1990081) q[2];
sx q[2];
rz(-1.4284632) q[2];
rz(-1.7510471) q[3];
sx q[3];
rz(-0.7775375) q[3];
sx q[3];
rz(-0.20974222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21935894) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(-1.6844164) q[0];
rz(2.4600785) q[1];
sx q[1];
rz(-0.69029713) q[1];
sx q[1];
rz(0.95726475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584314) q[0];
sx q[0];
rz(-0.51379075) q[0];
sx q[0];
rz(2.7130037) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8419015) q[2];
sx q[2];
rz(-0.77163358) q[2];
sx q[2];
rz(-2.6234805) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4217976) q[1];
sx q[1];
rz(-1.6198427) q[1];
sx q[1];
rz(-1.1625467) q[1];
rz(-pi) q[2];
rz(1.718347) q[3];
sx q[3];
rz(-0.06961623) q[3];
sx q[3];
rz(2.2152962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35931453) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(0.90685833) q[2];
rz(2.4685229) q[3];
sx q[3];
rz(-2.7570717) q[3];
sx q[3];
rz(-1.2955906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(0.54687706) q[0];
rz(-1.7710255) q[1];
sx q[1];
rz(-1.1616881) q[1];
sx q[1];
rz(-3.0414157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42368314) q[0];
sx q[0];
rz(-0.29742213) q[0];
sx q[0];
rz(-1.1017767) q[0];
rz(-pi) q[1];
rz(1.9630394) q[2];
sx q[2];
rz(-2.044495) q[2];
sx q[2];
rz(-2.365493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0292058) q[1];
sx q[1];
rz(-0.76923623) q[1];
sx q[1];
rz(-2.8112414) q[1];
rz(-pi) q[2];
rz(-1.1421575) q[3];
sx q[3];
rz(-1.1446307) q[3];
sx q[3];
rz(-1.7249191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.748041) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(2.098341) q[2];
rz(-2.4539963) q[3];
sx q[3];
rz(-2.6396773) q[3];
sx q[3];
rz(-0.25575328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461175) q[0];
sx q[0];
rz(-0.44545528) q[0];
sx q[0];
rz(-2.1050146) q[0];
rz(1.8702033) q[1];
sx q[1];
rz(-0.82415736) q[1];
sx q[1];
rz(-1.2870671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66726738) q[0];
sx q[0];
rz(-1.0762574) q[0];
sx q[0];
rz(0.70566341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.14187) q[2];
sx q[2];
rz(-2.1867036) q[2];
sx q[2];
rz(2.3162637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2894271) q[1];
sx q[1];
rz(-2.1680538) q[1];
sx q[1];
rz(0.35043033) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8356693) q[3];
sx q[3];
rz(-1.9960003) q[3];
sx q[3];
rz(-2.9028881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6048772) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(0.019901179) q[2];
rz(0.59602916) q[3];
sx q[3];
rz(-2.8607131) q[3];
sx q[3];
rz(1.9635487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6562281) q[0];
sx q[0];
rz(-1.2475659) q[0];
sx q[0];
rz(1.9670991) q[0];
rz(-0.30612048) q[1];
sx q[1];
rz(-2.0968292) q[1];
sx q[1];
rz(1.7721446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590483) q[0];
sx q[0];
rz(-1.567548) q[0];
sx q[0];
rz(-1.5643459) q[0];
rz(-2.9685176) q[2];
sx q[2];
rz(-2.1622116) q[2];
sx q[2];
rz(1.2966636) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8597653) q[1];
sx q[1];
rz(-1.9143811) q[1];
sx q[1];
rz(0.1600034) q[1];
x q[2];
rz(-1.9846652) q[3];
sx q[3];
rz(-1.0936495) q[3];
sx q[3];
rz(-2.8780495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29641637) q[2];
sx q[2];
rz(-2.3974819) q[2];
sx q[2];
rz(-2.3386492) q[2];
rz(-2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1882316) q[0];
sx q[0];
rz(-2.6022434) q[0];
sx q[0];
rz(-3.0292335) q[0];
rz(-0.27680963) q[1];
sx q[1];
rz(-1.9748297) q[1];
sx q[1];
rz(0.64754957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747276) q[0];
sx q[0];
rz(-1.234136) q[0];
sx q[0];
rz(-0.99322999) q[0];
rz(2.4729257) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(-0.47793717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6426443) q[1];
sx q[1];
rz(-2.0181839) q[1];
sx q[1];
rz(1.7018464) q[1];
rz(-1.3350491) q[3];
sx q[3];
rz(-1.7415348) q[3];
sx q[3];
rz(0.52908191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-0.15028149) q[2];
sx q[2];
rz(-1.3432937) q[2];
rz(0.51370931) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(0.59521365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419471) q[0];
sx q[0];
rz(-1.8193614) q[0];
sx q[0];
rz(2.7213726) q[0];
rz(1.3601903) q[1];
sx q[1];
rz(-1.9748634) q[1];
sx q[1];
rz(0.80418783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73690945) q[0];
sx q[0];
rz(-1.0361639) q[0];
sx q[0];
rz(-0.44628365) q[0];
rz(2.8560195) q[2];
sx q[2];
rz(-2.336506) q[2];
sx q[2];
rz(2.7325163) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11629352) q[1];
sx q[1];
rz(-2.5466047) q[1];
sx q[1];
rz(1.8354675) q[1];
rz(-pi) q[2];
rz(-3.0325417) q[3];
sx q[3];
rz(-1.4067459) q[3];
sx q[3];
rz(1.4993639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5899902) q[2];
sx q[2];
rz(-2.1253822) q[2];
sx q[2];
rz(1.972398) q[2];
rz(2.9679838) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(-1.9420697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46840295) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(-1.9386559) q[0];
rz(0.24083336) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(-2.208272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4543264) q[0];
sx q[0];
rz(-2.6250702) q[0];
sx q[0];
rz(2.3019426) q[0];
rz(-pi) q[1];
rz(-2.8779195) q[2];
sx q[2];
rz(-1.0142376) q[2];
sx q[2];
rz(-1.1900178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9874064) q[1];
sx q[1];
rz(-1.9498822) q[1];
sx q[1];
rz(-1.7861219) q[1];
rz(-pi) q[2];
rz(-2.1347789) q[3];
sx q[3];
rz(-1.6996592) q[3];
sx q[3];
rz(-3.0430678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11747083) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(0.13548279) q[2];
rz(-0.99902117) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(0.55618709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713772) q[0];
sx q[0];
rz(-1.1286292) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(1.0820092) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(-1.5622004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45856548) q[0];
sx q[0];
rz(-0.54046721) q[0];
sx q[0];
rz(1.6349273) q[0];
rz(0.57317971) q[2];
sx q[2];
rz(-0.64116353) q[2];
sx q[2];
rz(1.6784843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.434103) q[1];
sx q[1];
rz(-1.130778) q[1];
sx q[1];
rz(2.773519) q[1];
rz(-0.92023499) q[3];
sx q[3];
rz(-1.8585068) q[3];
sx q[3];
rz(0.77863151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18206) q[2];
sx q[2];
rz(-1.5675194) q[2];
sx q[2];
rz(0.62425295) q[2];
rz(-2.4652081) q[3];
sx q[3];
rz(-1.9703777) q[3];
sx q[3];
rz(-2.306365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312209) q[0];
sx q[0];
rz(-1.8932764) q[0];
sx q[0];
rz(-1.3488019) q[0];
rz(-2.8377332) q[1];
sx q[1];
rz(-1.5307348) q[1];
sx q[1];
rz(0.87711984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064909272) q[0];
sx q[0];
rz(-2.2936032) q[0];
sx q[0];
rz(-1.7641036) q[0];
rz(-pi) q[1];
rz(1.2563324) q[2];
sx q[2];
rz(-0.34032492) q[2];
sx q[2];
rz(-0.2257502) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7925229) q[1];
sx q[1];
rz(-1.9510837) q[1];
sx q[1];
rz(2.040634) q[1];
rz(-2.3043746) q[3];
sx q[3];
rz(-2.3851237) q[3];
sx q[3];
rz(-2.9314007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861686) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(1.0322734) q[2];
rz(-1.0885193) q[3];
sx q[3];
rz(-1.6531331) q[3];
sx q[3];
rz(0.77435875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9536204) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(-1.3183409) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(1.8404519) q[2];
sx q[2];
rz(-1.9741304) q[2];
sx q[2];
rz(-1.3000549) q[2];
rz(-1.5416038) q[3];
sx q[3];
rz(-0.43628449) q[3];
sx q[3];
rz(0.77602385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
