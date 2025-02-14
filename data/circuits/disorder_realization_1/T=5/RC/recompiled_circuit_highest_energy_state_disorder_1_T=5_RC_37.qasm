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
rz(-0.048592351) q[0];
sx q[0];
rz(-1.7162004) q[0];
sx q[0];
rz(-0.26560321) q[0];
rz(2.6939997) q[1];
sx q[1];
rz(-1.5040553) q[1];
sx q[1];
rz(0.12463364) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.671488) q[0];
sx q[0];
rz(-0.82297814) q[0];
sx q[0];
rz(2.1822071) q[0];
x q[1];
rz(1.9452478) q[2];
sx q[2];
rz(-1.6902897) q[2];
sx q[2];
rz(1.6993285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2217628) q[1];
sx q[1];
rz(-1.3897093) q[1];
sx q[1];
rz(-2.9403375) q[1];
rz(-pi) q[2];
rz(0.15480583) q[3];
sx q[3];
rz(-1.4034287) q[3];
sx q[3];
rz(-2.6449806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0521741) q[2];
sx q[2];
rz(-2.3201421) q[2];
sx q[2];
rz(0.87665147) q[2];
rz(-2.9540201) q[3];
sx q[3];
rz(-2.1552174) q[3];
sx q[3];
rz(-0.11876373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12980421) q[0];
sx q[0];
rz(-3.0846444) q[0];
sx q[0];
rz(2.7563128) q[0];
rz(2.724559) q[1];
sx q[1];
rz(-0.72154355) q[1];
sx q[1];
rz(1.0122274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4436662) q[0];
sx q[0];
rz(-1.3529643) q[0];
sx q[0];
rz(2.1429254) q[0];
rz(-pi) q[1];
rz(0.56372042) q[2];
sx q[2];
rz(-0.73838687) q[2];
sx q[2];
rz(0.2283048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9033563) q[1];
sx q[1];
rz(-0.84406772) q[1];
sx q[1];
rz(1.958117) q[1];
rz(-pi) q[2];
rz(-2.2506563) q[3];
sx q[3];
rz(-1.9417986) q[3];
sx q[3];
rz(1.551451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9932844) q[2];
sx q[2];
rz(-0.40996429) q[2];
sx q[2];
rz(-2.6170464) q[2];
rz(-2.2927393) q[3];
sx q[3];
rz(-0.69791228) q[3];
sx q[3];
rz(2.2822288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9808905) q[0];
sx q[0];
rz(-2.7301259) q[0];
sx q[0];
rz(-0.32660642) q[0];
rz(1.7713361) q[1];
sx q[1];
rz(-2.0854918) q[1];
sx q[1];
rz(-0.036570963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17509656) q[0];
sx q[0];
rz(-1.2155455) q[0];
sx q[0];
rz(0.47833832) q[0];
rz(-pi) q[1];
rz(2.6447949) q[2];
sx q[2];
rz(-2.3721176) q[2];
sx q[2];
rz(-2.0747607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72552437) q[1];
sx q[1];
rz(-1.6139133) q[1];
sx q[1];
rz(1.7579702) q[1];
rz(-pi) q[2];
rz(-2.1362334) q[3];
sx q[3];
rz(-0.32945368) q[3];
sx q[3];
rz(-2.0170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44564351) q[2];
sx q[2];
rz(-0.30569884) q[2];
sx q[2];
rz(1.1841527) q[2];
rz(0.46237692) q[3];
sx q[3];
rz(-2.0591044) q[3];
sx q[3];
rz(2.6760127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9584123) q[0];
sx q[0];
rz(-1.6540225) q[0];
sx q[0];
rz(2.2271449) q[0];
rz(0.86638266) q[1];
sx q[1];
rz(-1.2976846) q[1];
sx q[1];
rz(-0.31747174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7144517) q[0];
sx q[0];
rz(-1.1981816) q[0];
sx q[0];
rz(-1.6880549) q[0];
x q[1];
rz(-1.9201886) q[2];
sx q[2];
rz(-1.411806) q[2];
sx q[2];
rz(0.30245374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0054148) q[1];
sx q[1];
rz(-1.7247611) q[1];
sx q[1];
rz(1.2994159) q[1];
rz(-0.0013591493) q[3];
sx q[3];
rz(-1.1906644) q[3];
sx q[3];
rz(-1.330395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0471961) q[2];
sx q[2];
rz(-2.6698343) q[2];
sx q[2];
rz(2.6867234) q[2];
rz(1.9589849) q[3];
sx q[3];
rz(-2.9136361) q[3];
sx q[3];
rz(0.14127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32893786) q[0];
sx q[0];
rz(-1.8998572) q[0];
sx q[0];
rz(3.090233) q[0];
rz(-2.8142169) q[1];
sx q[1];
rz(-0.86762571) q[1];
sx q[1];
rz(0.059965722) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3695196) q[0];
sx q[0];
rz(-2.3926593) q[0];
sx q[0];
rz(-1.0706484) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9554638) q[2];
sx q[2];
rz(-0.95343381) q[2];
sx q[2];
rz(0.74194801) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0971036) q[1];
sx q[1];
rz(-1.8824235) q[1];
sx q[1];
rz(-2.1416041) q[1];
x q[2];
rz(-0.85598029) q[3];
sx q[3];
rz(-1.7188914) q[3];
sx q[3];
rz(0.87774897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7472234) q[2];
sx q[2];
rz(-2.2286434) q[2];
sx q[2];
rz(-2.8933914) q[2];
rz(0.40211755) q[3];
sx q[3];
rz(-0.44200236) q[3];
sx q[3];
rz(0.88002747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345271) q[0];
sx q[0];
rz(-1.2638673) q[0];
sx q[0];
rz(1.4136575) q[0];
rz(1.2029348) q[1];
sx q[1];
rz(-2.2212432) q[1];
sx q[1];
rz(3.038182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2884379) q[0];
sx q[0];
rz(-1.5010035) q[0];
sx q[0];
rz(-1.5091108) q[0];
rz(0.52210094) q[2];
sx q[2];
rz(-0.35873294) q[2];
sx q[2];
rz(-1.9821253) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6573665) q[1];
sx q[1];
rz(-1.5518589) q[1];
sx q[1];
rz(0.12469805) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9553071) q[3];
sx q[3];
rz(-1.8059192) q[3];
sx q[3];
rz(2.9540134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47205135) q[2];
sx q[2];
rz(-0.59110385) q[2];
sx q[2];
rz(-0.21311398) q[2];
rz(-1.537568) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(3.0311301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.9977426) q[0];
sx q[0];
rz(-2.282833) q[0];
sx q[0];
rz(2.6570038) q[0];
rz(-2.6735725) q[1];
sx q[1];
rz(-1.3222398) q[1];
sx q[1];
rz(-3.0048634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.384308) q[0];
sx q[0];
rz(-1.1955313) q[0];
sx q[0];
rz(2.882716) q[0];
x q[1];
rz(1.2263451) q[2];
sx q[2];
rz(-1.2023122) q[2];
sx q[2];
rz(2.0257575) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9715106) q[1];
sx q[1];
rz(-0.85935837) q[1];
sx q[1];
rz(2.9652639) q[1];
x q[2];
rz(2.8806512) q[3];
sx q[3];
rz(-0.43246864) q[3];
sx q[3];
rz(-2.3245426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4097269) q[2];
sx q[2];
rz(-0.12696433) q[2];
sx q[2];
rz(2.6597888) q[2];
rz(-1.9110154) q[3];
sx q[3];
rz(-1.1426208) q[3];
sx q[3];
rz(-0.13550152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9442673) q[0];
sx q[0];
rz(-2.8062286) q[0];
sx q[0];
rz(0.37099567) q[0];
rz(2.9839997) q[1];
sx q[1];
rz(-1.7540365) q[1];
sx q[1];
rz(-0.88034672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4828273) q[0];
sx q[0];
rz(-1.5012238) q[0];
sx q[0];
rz(2.0459904) q[0];
rz(-1.4029546) q[2];
sx q[2];
rz(-0.8693822) q[2];
sx q[2];
rz(-2.0244903) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0844565) q[1];
sx q[1];
rz(-1.5669654) q[1];
sx q[1];
rz(-1.4480897) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79788406) q[3];
sx q[3];
rz(-0.39643472) q[3];
sx q[3];
rz(2.1897922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2788435) q[2];
sx q[2];
rz(-2.4767488) q[2];
sx q[2];
rz(1.009602) q[2];
rz(-0.75262117) q[3];
sx q[3];
rz(-1.8372583) q[3];
sx q[3];
rz(0.93839222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143592) q[0];
sx q[0];
rz(-0.14150134) q[0];
sx q[0];
rz(-3.09521) q[0];
rz(1.4917689) q[1];
sx q[1];
rz(-1.7581419) q[1];
sx q[1];
rz(0.47328624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2641122) q[0];
sx q[0];
rz(-2.9213967) q[0];
sx q[0];
rz(-1.0408677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98433305) q[2];
sx q[2];
rz(-0.55078673) q[2];
sx q[2];
rz(0.065546602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1683732) q[1];
sx q[1];
rz(-2.1811317) q[1];
sx q[1];
rz(1.2471421) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1407309) q[3];
sx q[3];
rz(-1.5389082) q[3];
sx q[3];
rz(0.77939595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9778694) q[2];
sx q[2];
rz(-2.551584) q[2];
sx q[2];
rz(3.1150277) q[2];
rz(-0.50655347) q[3];
sx q[3];
rz(-2.3285464) q[3];
sx q[3];
rz(-0.51535362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.9764421) q[0];
sx q[0];
rz(-0.9557752) q[0];
sx q[0];
rz(2.9859848) q[0];
rz(-2.7071629) q[1];
sx q[1];
rz(-2.5948718) q[1];
sx q[1];
rz(0.2880407) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018943448) q[0];
sx q[0];
rz(-0.83479202) q[0];
sx q[0];
rz(-2.3998866) q[0];
x q[1];
rz(1.4819809) q[2];
sx q[2];
rz(-0.15379158) q[2];
sx q[2];
rz(1.8078992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0295683) q[1];
sx q[1];
rz(-2.3178737) q[1];
sx q[1];
rz(2.1805641) q[1];
rz(-pi) q[2];
rz(-2.1197917) q[3];
sx q[3];
rz(-2.4394301) q[3];
sx q[3];
rz(0.3265061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7274999) q[2];
sx q[2];
rz(-0.45656559) q[2];
sx q[2];
rz(-3.0957733) q[2];
rz(0.11848655) q[3];
sx q[3];
rz(-0.89209569) q[3];
sx q[3];
rz(-2.4011325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988149) q[0];
sx q[0];
rz(-1.6294263) q[0];
sx q[0];
rz(-1.908041) q[0];
rz(-1.9998101) q[1];
sx q[1];
rz(-0.70527609) q[1];
sx q[1];
rz(-1.3377778) q[1];
rz(-2.4547443) q[2];
sx q[2];
rz(-2.3301132) q[2];
sx q[2];
rz(-2.2215497) q[2];
rz(-1.508504) q[3];
sx q[3];
rz(-1.6332165) q[3];
sx q[3];
rz(-1.6467057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
