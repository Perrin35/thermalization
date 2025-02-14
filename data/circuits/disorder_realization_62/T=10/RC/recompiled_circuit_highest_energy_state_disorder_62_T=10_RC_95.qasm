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
rz(-2.488774) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(0.62243661) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1718423) q[0];
sx q[0];
rz(-1.4542219) q[0];
sx q[0];
rz(1.2571093) q[0];
rz(-pi) q[1];
rz(1.2369701) q[2];
sx q[2];
rz(-0.99572748) q[2];
sx q[2];
rz(1.0525557) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37934946) q[1];
sx q[1];
rz(-2.6062638) q[1];
sx q[1];
rz(-1.204727) q[1];
x q[2];
rz(-2.1795737) q[3];
sx q[3];
rz(-1.7005657) q[3];
sx q[3];
rz(-2.6994841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6301253) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-2.5959065) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-2.4617709) q[3];
sx q[3];
rz(-0.95664501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086455258) q[0];
sx q[0];
rz(-2.4392023) q[0];
sx q[0];
rz(-2.0461244) q[0];
rz(-0.99758863) q[1];
sx q[1];
rz(-1.2350524) q[1];
sx q[1];
rz(1.197061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8923949) q[0];
sx q[0];
rz(-0.64606842) q[0];
sx q[0];
rz(1.9996253) q[0];
rz(-pi) q[1];
rz(-0.95112339) q[2];
sx q[2];
rz(-0.89604267) q[2];
sx q[2];
rz(-2.5017966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61651308) q[1];
sx q[1];
rz(-2.5933806) q[1];
sx q[1];
rz(2.1384778) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6359699) q[3];
sx q[3];
rz(-2.7693336) q[3];
sx q[3];
rz(0.79141189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41584388) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.7572629) q[2];
rz(-1.7978801) q[3];
sx q[3];
rz(-1.5680771) q[3];
sx q[3];
rz(1.869092) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5196853) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(-2.7144879) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.4991263) q[1];
sx q[1];
rz(-2.5274091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489036) q[0];
sx q[0];
rz(-1.2430191) q[0];
sx q[0];
rz(1.8564188) q[0];
rz(-pi) q[1];
rz(0.74684192) q[2];
sx q[2];
rz(-1.421442) q[2];
sx q[2];
rz(-2.9666882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79518049) q[1];
sx q[1];
rz(-1.507057) q[1];
sx q[1];
rz(1.6883114) q[1];
rz(-pi) q[2];
rz(-2.6355686) q[3];
sx q[3];
rz(-1.783924) q[3];
sx q[3];
rz(-2.0863078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68982879) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(-2.8705987) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52151742) q[0];
sx q[0];
rz(-1.2780315) q[0];
sx q[0];
rz(-1.0241246) q[0];
x q[1];
rz(-0.35778862) q[2];
sx q[2];
rz(-2.3364106) q[2];
sx q[2];
rz(1.5186269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4252538) q[1];
sx q[1];
rz(-2.0074866) q[1];
sx q[1];
rz(1.727987) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3621394) q[3];
sx q[3];
rz(-1.28445) q[3];
sx q[3];
rz(-1.7735667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0858687) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(1.9191939) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(-0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5772575) q[0];
sx q[0];
rz(-0.99522796) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(1.0352146) q[1];
sx q[1];
rz(-1.2151006) q[1];
sx q[1];
rz(1.530102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852762) q[0];
sx q[0];
rz(-1.0438598) q[0];
sx q[0];
rz(-2.2599392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3283073) q[2];
sx q[2];
rz(-1.7264778) q[2];
sx q[2];
rz(-0.98634431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4799696) q[1];
sx q[1];
rz(-1.1692059) q[1];
sx q[1];
rz(1.5127403) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6282998) q[3];
sx q[3];
rz(-1.1604476) q[3];
sx q[3];
rz(-1.8156605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(-0.063892603) q[2];
rz(-0.70819267) q[3];
sx q[3];
rz(-1.4539366) q[3];
sx q[3];
rz(-1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42234364) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(-1.0759906) q[0];
rz(-0.05323449) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(2.7630973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196825) q[0];
sx q[0];
rz(-1.440597) q[0];
sx q[0];
rz(0.2567592) q[0];
rz(-pi) q[1];
rz(-3.1000146) q[2];
sx q[2];
rz(-1.4828621) q[2];
sx q[2];
rz(-1.8454332) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4617745) q[1];
sx q[1];
rz(-1.8289806) q[1];
sx q[1];
rz(-1.3854909) q[1];
rz(-pi) q[2];
rz(0.58121292) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(0.17457822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1899015) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(-1.0339197) q[2];
rz(-0.90384358) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(-3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251625) q[0];
sx q[0];
rz(-1.4465541) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(-0.12475363) q[1];
sx q[1];
rz(-1.1589103) q[1];
sx q[1];
rz(-0.4932901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.497555) q[0];
sx q[0];
rz(-1.8480267) q[0];
sx q[0];
rz(-1.9009186) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5458844) q[2];
sx q[2];
rz(-1.7988211) q[2];
sx q[2];
rz(1.4268503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64831454) q[1];
sx q[1];
rz(-1.1469134) q[1];
sx q[1];
rz(-2.5878169) q[1];
rz(-pi) q[2];
rz(0.36549632) q[3];
sx q[3];
rz(-0.29416829) q[3];
sx q[3];
rz(-0.95137596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7571681) q[2];
sx q[2];
rz(-1.9809456) q[2];
sx q[2];
rz(1.9942795) q[2];
rz(0.82289639) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(-2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-0.4185032) q[0];
rz(3.0283527) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(2.0638827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5302056) q[0];
sx q[0];
rz(-0.28908397) q[0];
sx q[0];
rz(1.6340088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8968717) q[2];
sx q[2];
rz(-0.80279826) q[2];
sx q[2];
rz(1.7434415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7188765) q[1];
sx q[1];
rz(-0.36564974) q[1];
sx q[1];
rz(-1.7083733) q[1];
x q[2];
rz(-0.65656482) q[3];
sx q[3];
rz(-1.5235008) q[3];
sx q[3];
rz(1.8790203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30617994) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(-1.6861247) q[2];
rz(0.16737394) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48285943) q[0];
sx q[0];
rz(-1.3128244) q[0];
sx q[0];
rz(0.40153781) q[0];
rz(2.7032779) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(1.6848791) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4171203) q[0];
sx q[0];
rz(-0.03532413) q[0];
sx q[0];
rz(1.2958584) q[0];
x q[1];
rz(-2.8598665) q[2];
sx q[2];
rz(-0.67517074) q[2];
sx q[2];
rz(-0.074568579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40565421) q[1];
sx q[1];
rz(-0.88684139) q[1];
sx q[1];
rz(2.4962462) q[1];
x q[2];
rz(-2.4026186) q[3];
sx q[3];
rz(-2.4904076) q[3];
sx q[3];
rz(0.84924752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38176408) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(-2.4299202) q[2];
rz(0.034218637) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(-1.155352) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3592247) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(2.0955739) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(-2.9050713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6470784) q[0];
sx q[0];
rz(-0.064726742) q[0];
sx q[0];
rz(-1.2184124) q[0];
rz(-2.3600618) q[2];
sx q[2];
rz(-2.3252914) q[2];
sx q[2];
rz(-2.3280509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1139086) q[1];
sx q[1];
rz(-2.5792173) q[1];
sx q[1];
rz(-0.67090583) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-0.5566842) q[2];
sx q[2];
rz(0.64209783) q[2];
rz(0.56946483) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6713329) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.9500465) q[1];
sx q[1];
rz(-1.5622495) q[1];
sx q[1];
rz(-1.539485) q[1];
rz(-1.446522) q[2];
sx q[2];
rz(-2.3588603) q[2];
sx q[2];
rz(-1.0609577) q[2];
rz(-2.8068918) q[3];
sx q[3];
rz(-1.1664433) q[3];
sx q[3];
rz(0.56952624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
