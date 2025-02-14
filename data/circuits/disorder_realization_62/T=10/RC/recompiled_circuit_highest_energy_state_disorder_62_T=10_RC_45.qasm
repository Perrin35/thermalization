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
rz(-0.57617968) q[0];
sx q[0];
rz(-1.4644858) q[0];
sx q[0];
rz(-0.65281868) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(-2.519156) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718423) q[0];
sx q[0];
rz(-1.4542219) q[0];
sx q[0];
rz(1.8844834) q[0];
x q[1];
rz(-1.2369701) q[2];
sx q[2];
rz(-0.99572748) q[2];
sx q[2];
rz(-1.0525557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7622432) q[1];
sx q[1];
rz(-2.6062638) q[1];
sx q[1];
rz(-1.204727) q[1];
rz(-0.15776026) q[3];
sx q[3];
rz(-0.96786849) q[3];
sx q[3];
rz(-1.2186528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6301253) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(0.54568616) q[2];
rz(-0.68583471) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(-2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551374) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(-2.0461244) q[0];
rz(0.99758863) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.9445317) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491978) q[0];
sx q[0];
rz(-0.64606842) q[0];
sx q[0];
rz(-1.1419673) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77669528) q[2];
sx q[2];
rz(-2.0414845) q[2];
sx q[2];
rz(-1.7913417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8833739) q[1];
sx q[1];
rz(-1.1158459) q[1];
sx q[1];
rz(-2.8243498) q[1];
rz(-pi) q[2];
rz(-1.7577111) q[3];
sx q[3];
rz(-1.2469562) q[3];
sx q[3];
rz(1.8139386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41584388) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.3843298) q[2];
rz(-1.3437126) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(-1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5196853) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(0.4271048) q[0];
rz(-1.2318132) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(-2.5274091) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489036) q[0];
sx q[0];
rz(-1.8985735) q[0];
sx q[0];
rz(-1.8564188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3947507) q[2];
sx q[2];
rz(-1.421442) q[2];
sx q[2];
rz(-0.17490444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28084785) q[1];
sx q[1];
rz(-3.0079746) q[1];
sx q[1];
rz(1.0722776) q[1];
rz(-0.50602405) q[3];
sx q[3];
rz(-1.783924) q[3];
sx q[3];
rz(-1.0552849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4517639) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(1.8577925) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2778306) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(-2.7226287) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(3.0640501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5351535) q[0];
sx q[0];
rz(-0.61302671) q[0];
sx q[0];
rz(1.0453348) q[0];
rz(-pi) q[1];
x q[1];
rz(2.783804) q[2];
sx q[2];
rz(-2.3364106) q[2];
sx q[2];
rz(-1.6229657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0748914) q[1];
sx q[1];
rz(-2.6791926) q[1];
sx q[1];
rz(-2.8179864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7794533) q[3];
sx q[3];
rz(-1.8571427) q[3];
sx q[3];
rz(1.3680259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.055723995) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(-0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(1.0352146) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.530102) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28883067) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(-2.2599392) q[0];
x q[1];
rz(0.45299977) q[2];
sx q[2];
rz(-0.36213798) q[2];
sx q[2];
rz(1.0116133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6275806) q[1];
sx q[1];
rz(-0.40553933) q[1];
sx q[1];
rz(-0.13579129) q[1];
rz(0.13134457) q[3];
sx q[3];
rz(-2.7274611) q[3];
sx q[3];
rz(-1.9589748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3349541) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(-3.0777001) q[2];
rz(-0.70819267) q[3];
sx q[3];
rz(-1.4539366) q[3];
sx q[3];
rz(1.8074869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42234364) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(-1.0759906) q[0];
rz(-0.05323449) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(2.7630973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16679487) q[0];
sx q[0];
rz(-2.85436) q[0];
sx q[0];
rz(2.6655281) q[0];
x q[1];
rz(-3.1000146) q[2];
sx q[2];
rz(-1.4828621) q[2];
sx q[2];
rz(-1.8454332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9847451) q[1];
sx q[1];
rz(-1.7498921) q[1];
sx q[1];
rz(-2.8791134) q[1];
x q[2];
rz(0.58121292) q[3];
sx q[3];
rz(-2.2848115) q[3];
sx q[3];
rz(-0.17457822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1899015) q[2];
sx q[2];
rz(-1.6804164) q[2];
sx q[2];
rz(-2.107673) q[2];
rz(-2.2377491) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(-3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164301) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(-3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(-0.4932901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64403764) q[0];
sx q[0];
rz(-1.8480267) q[0];
sx q[0];
rz(1.2406741) q[0];
rz(-1.5957082) q[2];
sx q[2];
rz(-1.3427715) q[2];
sx q[2];
rz(-1.4268503) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3354937) q[1];
sx q[1];
rz(-0.68365288) q[1];
sx q[1];
rz(-2.4324576) q[1];
x q[2];
rz(-0.27573891) q[3];
sx q[3];
rz(-1.4669802) q[3];
sx q[3];
rz(2.8732515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3844246) q[2];
sx q[2];
rz(-1.9809456) q[2];
sx q[2];
rz(-1.9942795) q[2];
rz(-0.82289639) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(-0.66250044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61854521) q[0];
sx q[0];
rz(-1.2724027) q[0];
sx q[0];
rz(2.7230895) q[0];
rz(0.11323994) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(-2.0638827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5302056) q[0];
sx q[0];
rz(-2.8525087) q[0];
sx q[0];
rz(1.6340088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.244721) q[2];
sx q[2];
rz(-2.3387944) q[2];
sx q[2];
rz(1.3981512) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.019494836) q[1];
sx q[1];
rz(-1.6198525) q[1];
sx q[1];
rz(8*pi/13) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.077385888) q[3];
sx q[3];
rz(-2.4835781) q[3];
sx q[3];
rz(-0.24694996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30617994) q[2];
sx q[2];
rz(-0.70859185) q[2];
sx q[2];
rz(-1.4554679) q[2];
rz(-0.16737394) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(-1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48285943) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(0.40153781) q[0];
rz(-0.43831476) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.4567136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4171203) q[0];
sx q[0];
rz(-3.1062685) q[0];
sx q[0];
rz(-1.2958584) q[0];
x q[1];
rz(-2.8598665) q[2];
sx q[2];
rz(-2.4664219) q[2];
sx q[2];
rz(0.074568579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46692013) q[1];
sx q[1];
rz(-0.90306696) q[1];
sx q[1];
rz(0.93514644) q[1];
x q[2];
rz(2.4026186) q[3];
sx q[3];
rz(-0.65118507) q[3];
sx q[3];
rz(0.84924752) q[3];
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
rz(-0.70845571) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.782368) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(0.79886287) q[0];
rz(-2.0955739) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(2.9050713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49451429) q[0];
sx q[0];
rz(-0.064726742) q[0];
sx q[0];
rz(-1.2184124) q[0];
rz(-0.78153083) q[2];
sx q[2];
rz(-2.3252914) q[2];
sx q[2];
rz(-0.81354173) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1139086) q[1];
sx q[1];
rz(-2.5792173) q[1];
sx q[1];
rz(2.4706868) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71098401) q[3];
sx q[3];
rz(-1.4812638) q[3];
sx q[3];
rz(0.79345616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.11655) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(2.4994948) q[2];
rz(-0.56946483) q[3];
sx q[3];
rz(-0.48782188) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.1915462) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(-1.446522) q[2];
sx q[2];
rz(-2.3588603) q[2];
sx q[2];
rz(-1.0609577) q[2];
rz(2.8068918) q[3];
sx q[3];
rz(-1.9751493) q[3];
sx q[3];
rz(-2.5720664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
