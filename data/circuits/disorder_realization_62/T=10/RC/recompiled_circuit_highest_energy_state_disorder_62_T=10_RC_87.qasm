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
rz(2.565413) q[0];
sx q[0];
rz(-1.6771069) q[0];
sx q[0];
rz(0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(5.6607487) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718423) q[0];
sx q[0];
rz(-1.4542219) q[0];
sx q[0];
rz(1.8844834) q[0];
rz(-2.6735252) q[2];
sx q[2];
rz(-0.65535802) q[2];
sx q[2];
rz(1.6201902) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6316526) q[1];
sx q[1];
rz(-1.7544244) q[1];
sx q[1];
rz(1.0650403) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1795737) q[3];
sx q[3];
rz(-1.441027) q[3];
sx q[3];
rz(2.6994841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5114674) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-0.54568616) q[2];
rz(0.68583471) q[3];
sx q[3];
rz(-2.4617709) q[3];
sx q[3];
rz(0.95664501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086455258) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(-2.0461244) q[0];
rz(0.99758863) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.9445317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491978) q[0];
sx q[0];
rz(-0.64606842) q[0];
sx q[0];
rz(-1.1419673) q[0];
rz(-0.77669528) q[2];
sx q[2];
rz(-1.1001081) q[2];
sx q[2];
rz(1.7913417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8833739) q[1];
sx q[1];
rz(-2.0257468) q[1];
sx q[1];
rz(-0.31724288) q[1];
x q[2];
rz(0.50562276) q[3];
sx q[3];
rz(-0.37225906) q[3];
sx q[3];
rz(2.3501808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41584388) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.3843298) q[2];
rz(1.7978801) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(-1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62190732) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(-2.7144879) q[0];
rz(1.9097795) q[1];
sx q[1];
rz(-1.4991263) q[1];
sx q[1];
rz(2.5274091) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012497525) q[0];
sx q[0];
rz(-1.84082) q[0];
sx q[0];
rz(-0.34059033) q[0];
rz(-pi) q[1];
x q[1];
rz(1.368548) q[2];
sx q[2];
rz(-0.83421153) q[2];
sx q[2];
rz(-1.5327765) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28084785) q[1];
sx q[1];
rz(-3.0079746) q[1];
sx q[1];
rz(-2.0693151) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3282458) q[3];
sx q[3];
rz(-1.0772675) q[3];
sx q[3];
rz(-0.39883116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68982879) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(2.8705987) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2778306) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(0.41896391) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(0.077542543) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5351535) q[0];
sx q[0];
rz(-0.61302671) q[0];
sx q[0];
rz(-1.0453348) q[0];
rz(-pi) q[1];
rz(-1.9201905) q[2];
sx q[2];
rz(-2.3121693) q[2];
sx q[2];
rz(-1.0238613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71633881) q[1];
sx q[1];
rz(-2.0074866) q[1];
sx q[1];
rz(1.727987) q[1];
rz(-2.8492498) q[3];
sx q[3];
rz(-1.3707531) q[3];
sx q[3];
rz(0.26250473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.055723995) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(1.2223988) q[2];
rz(-2.870765) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643352) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.3539535) q[0];
rz(-1.0352146) q[1];
sx q[1];
rz(-1.2151006) q[1];
sx q[1];
rz(1.6114906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852762) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(-0.88165347) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6885929) q[2];
sx q[2];
rz(-0.36213798) q[2];
sx q[2];
rz(2.1299794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93188796) q[1];
sx q[1];
rz(-1.5173638) q[1];
sx q[1];
rz(0.40219743) q[1];
x q[2];
rz(1.5132929) q[3];
sx q[3];
rz(-1.981145) q[3];
sx q[3];
rz(-1.8156605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(-3.0777001) q[2];
rz(2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719249) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(-1.0759906) q[0];
rz(-0.05323449) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(-0.3784953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16679487) q[0];
sx q[0];
rz(-2.85436) q[0];
sx q[0];
rz(-0.47606456) q[0];
x q[1];
rz(0.041578023) q[2];
sx q[2];
rz(-1.4828621) q[2];
sx q[2];
rz(1.2961594) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15684756) q[1];
sx q[1];
rz(-1.7498921) q[1];
sx q[1];
rz(-2.8791134) q[1];
rz(-2.3742484) q[3];
sx q[3];
rz(-1.9986614) q[3];
sx q[3];
rz(-1.3391173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1899015) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(1.0339197) q[2];
rz(0.90384358) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(-3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164301) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(-2.548581) q[0];
rz(3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(-2.6483026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.497555) q[0];
sx q[0];
rz(-1.293566) q[0];
sx q[0];
rz(1.9009186) q[0];
rz(0.10693018) q[2];
sx q[2];
rz(-2.9122346) q[2];
sx q[2];
rz(1.3170674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3354937) q[1];
sx q[1];
rz(-0.68365288) q[1];
sx q[1];
rz(0.70913507) q[1];
x q[2];
rz(-0.36549632) q[3];
sx q[3];
rz(-2.8474244) q[3];
sx q[3];
rz(2.1902167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7571681) q[2];
sx q[2];
rz(-1.9809456) q[2];
sx q[2];
rz(-1.9942795) q[2];
rz(-2.3186963) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(0.66250044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(2.7230895) q[0];
rz(3.0283527) q[1];
sx q[1];
rz(-1.4694045) q[1];
sx q[1];
rz(-2.0638827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200054) q[0];
sx q[0];
rz(-1.5888056) q[0];
sx q[0];
rz(1.8593345) q[0];
x q[1];
rz(2.3465326) q[2];
sx q[2];
rz(-1.3382947) q[2];
sx q[2];
rz(2.7382221) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.019494836) q[1];
sx q[1];
rz(-1.5217402) q[1];
sx q[1];
rz(-8*pi/13) q[1];
rz(-pi) q[2];
rz(-0.65656482) q[3];
sx q[3];
rz(-1.5235008) q[3];
sx q[3];
rz(-1.2625723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8354127) q[2];
sx q[2];
rz(-0.70859185) q[2];
sx q[2];
rz(-1.4554679) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(-2.0696056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48285943) q[0];
sx q[0];
rz(-1.3128244) q[0];
sx q[0];
rz(0.40153781) q[0];
rz(0.43831476) q[1];
sx q[1];
rz(-0.79151789) q[1];
sx q[1];
rz(-1.4567136) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4171203) q[0];
sx q[0];
rz(-0.03532413) q[0];
sx q[0];
rz(-1.8457343) q[0];
rz(-1.3517604) q[2];
sx q[2];
rz(-2.214785) q[2];
sx q[2];
rz(-0.42962675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5324133) q[1];
sx q[1];
rz(-1.0858156) q[1];
sx q[1];
rz(0.77528016) q[1];
rz(-0.73897408) q[3];
sx q[3];
rz(-0.65118507) q[3];
sx q[3];
rz(-2.2923451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7598286) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(-0.71167243) q[2];
rz(0.034218637) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(1.9862407) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782368) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(2.0955739) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(2.9050713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49451429) q[0];
sx q[0];
rz(-0.064726742) q[0];
sx q[0];
rz(-1.9231803) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2138519) q[2];
sx q[2];
rz(-2.1143713) q[2];
sx q[2];
rz(1.7801931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7812209) q[1];
sx q[1];
rz(-1.1399593) q[1];
sx q[1];
rz(-1.9442417) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71098401) q[3];
sx q[3];
rz(-1.6603289) q[3];
sx q[3];
rz(-0.79345616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(2.4994948) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.9500465) q[1];
sx q[1];
rz(-1.5622495) q[1];
sx q[1];
rz(-1.539485) q[1];
rz(-2.3496579) q[2];
sx q[2];
rz(-1.4832693) q[2];
sx q[2];
rz(-2.5434189) q[2];
rz(-0.33470086) q[3];
sx q[3];
rz(-1.9751493) q[3];
sx q[3];
rz(-2.5720664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
