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
rz(-0.0039761877) q[0];
sx q[0];
rz(-3.0206326) q[0];
rz(0.51789415) q[1];
sx q[1];
rz(-0.33101141) q[1];
sx q[1];
rz(-1.1878699) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8594583) q[0];
sx q[0];
rz(-1.6573424) q[0];
sx q[0];
rz(0.29973932) q[0];
rz(1.3224363) q[2];
sx q[2];
rz(-2.6018086) q[2];
sx q[2];
rz(2.7960461) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29562274) q[1];
sx q[1];
rz(-0.69087183) q[1];
sx q[1];
rz(0.9102896) q[1];
rz(-1.6191606) q[3];
sx q[3];
rz(-0.93975583) q[3];
sx q[3];
rz(3.0121355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0594222) q[2];
sx q[2];
rz(-0.94258451) q[2];
sx q[2];
rz(-1.7131294) q[2];
rz(-1.7510471) q[3];
sx q[3];
rz(-0.7775375) q[3];
sx q[3];
rz(-0.20974222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9222337) q[0];
sx q[0];
rz(-1.5606422) q[0];
sx q[0];
rz(-1.4571762) q[0];
rz(-2.4600785) q[1];
sx q[1];
rz(-2.4512955) q[1];
sx q[1];
rz(0.95726475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5752561) q[0];
sx q[0];
rz(-2.0342376) q[0];
sx q[0];
rz(-1.3404247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.886734) q[2];
sx q[2];
rz(-2.3074842) q[2];
sx q[2];
rz(-2.2535549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87220472) q[1];
sx q[1];
rz(-1.1630668) q[1];
sx q[1];
rz(0.053430037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5019341) q[3];
sx q[3];
rz(-1.5810229) q[3];
sx q[3];
rz(2.6442912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35931453) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(-2.2347343) q[2];
rz(0.67306972) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(-1.2955906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(2.5947156) q[0];
rz(-1.3705672) q[1];
sx q[1];
rz(-1.9799045) q[1];
sx q[1];
rz(0.10017698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4456383) q[0];
sx q[0];
rz(-1.7036519) q[0];
sx q[0];
rz(-1.3039051) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1785533) q[2];
sx q[2];
rz(-2.044495) q[2];
sx q[2];
rz(-2.365493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11238687) q[1];
sx q[1];
rz(-0.76923623) q[1];
sx q[1];
rz(2.8112414) q[1];
x q[2];
rz(-2.4002714) q[3];
sx q[3];
rz(-2.5466515) q[3];
sx q[3];
rz(-2.5603608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.748041) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(-2.098341) q[2];
rz(-2.4539963) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(0.25575328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8461175) q[0];
sx q[0];
rz(-0.44545528) q[0];
sx q[0];
rz(2.1050146) q[0];
rz(1.2713894) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(1.8545256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51924473) q[0];
sx q[0];
rz(-2.17832) q[0];
sx q[0];
rz(-2.1871479) q[0];
rz(-pi) q[1];
rz(0.66133581) q[2];
sx q[2];
rz(-1.2244389) q[2];
sx q[2];
rz(2.6544184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85216552) q[1];
sx q[1];
rz(-2.1680538) q[1];
sx q[1];
rz(0.35043033) q[1];
rz(-pi) q[2];
rz(2.1577094) q[3];
sx q[3];
rz(-2.623284) q[3];
sx q[3];
rz(-0.89215088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6048772) q[2];
sx q[2];
rz(-0.36467364) q[2];
sx q[2];
rz(0.019901179) q[2];
rz(2.5455635) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(1.9635487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-0.48536456) q[0];
sx q[0];
rz(-1.2475659) q[0];
sx q[0];
rz(-1.9670991) q[0];
rz(-2.8354722) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(1.7721446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4547538) q[0];
sx q[0];
rz(-3.1343705) q[0];
sx q[0];
rz(-1.104284) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9685176) q[2];
sx q[2];
rz(-0.97938108) q[2];
sx q[2];
rz(-1.8449291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9766531) q[1];
sx q[1];
rz(-0.37766981) q[1];
sx q[1];
rz(1.151848) q[1];
rz(-pi) q[2];
rz(2.4804875) q[3];
sx q[3];
rz(-0.62097681) q[3];
sx q[3];
rz(-1.0266765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8451763) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(2.3386492) q[2];
rz(-1.1154277) q[3];
sx q[3];
rz(-2.4468456) q[3];
sx q[3];
rz(1.7723134) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95336103) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(0.11235919) q[0];
rz(-0.27680963) q[1];
sx q[1];
rz(-1.1667629) q[1];
sx q[1];
rz(-0.64754957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43418068) q[0];
sx q[0];
rz(-0.65874225) q[0];
sx q[0];
rz(1.0007658) q[0];
x q[1];
rz(0.66866691) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(-2.6636555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9383835) q[1];
sx q[1];
rz(-2.6766549) q[1];
sx q[1];
rz(0.26588456) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3350491) q[3];
sx q[3];
rz(-1.4000579) q[3];
sx q[3];
rz(-2.6125107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67927805) q[2];
sx q[2];
rz(-0.15028149) q[2];
sx q[2];
rz(-1.3432937) q[2];
rz(0.51370931) q[3];
sx q[3];
rz(-1.4880344) q[3];
sx q[3];
rz(2.546379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419471) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(0.42022002) q[0];
rz(1.3601903) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(-0.80418783) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4046832) q[0];
sx q[0];
rz(-1.0361639) q[0];
sx q[0];
rz(2.695309) q[0];
rz(2.3571722) q[2];
sx q[2];
rz(-1.7752991) q[2];
sx q[2];
rz(0.96098268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7090164) q[1];
sx q[1];
rz(-2.1423999) q[1];
sx q[1];
rz(0.17523058) q[1];
rz(-pi) q[2];
rz(3.0325417) q[3];
sx q[3];
rz(-1.4067459) q[3];
sx q[3];
rz(-1.4993639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5516025) q[2];
sx q[2];
rz(-1.0162105) q[2];
sx q[2];
rz(-1.1691947) q[2];
rz(0.17360887) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(1.9420697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(2.6731897) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(-1.9386559) q[0];
rz(2.9007593) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(2.208272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.653395) q[0];
sx q[0];
rz(-1.9472593) q[0];
sx q[0];
rz(-2.7791136) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9674932) q[2];
sx q[2];
rz(-2.5317414) q[2];
sx q[2];
rz(1.4791453) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62040239) q[1];
sx q[1];
rz(-0.4333843) q[1];
sx q[1];
rz(-0.49232011) q[1];
rz(-pi) q[2];
rz(-1.0068137) q[3];
sx q[3];
rz(-1.6996592) q[3];
sx q[3];
rz(3.0430678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0241218) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(3.0061099) q[2];
rz(0.99902117) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67021543) q[0];
sx q[0];
rz(-2.0129634) q[0];
sx q[0];
rz(-1.7673329) q[0];
rz(1.0820092) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(-1.5622004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53331965) q[0];
sx q[0];
rz(-1.0315622) q[0];
sx q[0];
rz(0.03843741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5815046) q[2];
sx q[2];
rz(-1.9011407) q[2];
sx q[2];
rz(0.3696839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.434103) q[1];
sx q[1];
rz(-2.0108147) q[1];
sx q[1];
rz(-2.773519) q[1];
rz(-pi) q[2];
rz(0.92023499) q[3];
sx q[3];
rz(-1.2830858) q[3];
sx q[3];
rz(-2.3629611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18206) q[2];
sx q[2];
rz(-1.5740732) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4312209) q[0];
sx q[0];
rz(-1.8932764) q[0];
sx q[0];
rz(1.3488019) q[0];
rz(-2.8377332) q[1];
sx q[1];
rz(-1.5307348) q[1];
sx q[1];
rz(0.87711984) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9187885) q[0];
sx q[0];
rz(-2.3979146) q[0];
sx q[0];
rz(-0.2144465) q[0];
rz(-pi) q[1];
rz(-1.8956081) q[2];
sx q[2];
rz(-1.6742252) q[2];
sx q[2];
rz(-1.6425486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7925229) q[1];
sx q[1];
rz(-1.9510837) q[1];
sx q[1];
rz(-2.040634) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3043746) q[3];
sx q[3];
rz(-0.756469) q[3];
sx q[3];
rz(-2.9314007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0861686) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(1.0322734) q[2];
rz(-1.0885193) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(-0.77435875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9536204) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(1.3183409) q[1];
sx q[1];
rz(-1.4069469) q[1];
sx q[1];
rz(-2.5896752) q[1];
rz(1.3011408) q[2];
sx q[2];
rz(-1.1674623) q[2];
sx q[2];
rz(1.8415378) q[2];
rz(-0.013608215) q[3];
sx q[3];
rz(-1.1347105) q[3];
sx q[3];
rz(-2.3977765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
