OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(3.1443449) q[0];
sx q[0];
rz(10.436463) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(-1.8399532) q[1];
sx q[1];
rz(-0.17170061) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8780706) q[0];
sx q[0];
rz(-1.6524914) q[0];
sx q[0];
rz(-0.6058713) q[0];
rz(-pi) q[1];
rz(0.5131716) q[2];
sx q[2];
rz(-0.22226873) q[2];
sx q[2];
rz(-1.5813511) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2312647) q[1];
sx q[1];
rz(-0.99580169) q[1];
sx q[1];
rz(1.0564338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0461476) q[3];
sx q[3];
rz(-0.10766115) q[3];
sx q[3];
rz(-0.3357418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95919886) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(-1.4200776) q[2];
rz(-1.1848263) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0576393) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-0.49348304) q[0];
rz(0.77589846) q[1];
sx q[1];
rz(-2.6319365) q[1];
sx q[1];
rz(2.6367771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579266) q[0];
sx q[0];
rz(-0.89910451) q[0];
sx q[0];
rz(-2.2855285) q[0];
rz(-pi) q[1];
rz(2.6146982) q[2];
sx q[2];
rz(-1.3119446) q[2];
sx q[2];
rz(-2.7833454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.069020443) q[1];
sx q[1];
rz(-2.1694599) q[1];
sx q[1];
rz(2.0607949) q[1];
rz(-pi) q[2];
rz(-2.9166031) q[3];
sx q[3];
rz(-2.8691022) q[3];
sx q[3];
rz(2.837473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8043171) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(2.5110631) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(1.4001018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0640963) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(-0.65573829) q[0];
rz(2.8302622) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(3.0701367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27256672) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-0.053215543) q[0];
rz(-pi) q[1];
rz(-1.9391483) q[2];
sx q[2];
rz(-1.2542274) q[2];
sx q[2];
rz(-2.9376415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7332257) q[1];
sx q[1];
rz(-1.0828312) q[1];
sx q[1];
rz(-0.074175149) q[1];
rz(-pi) q[2];
rz(-0.53160588) q[3];
sx q[3];
rz(-1.387145) q[3];
sx q[3];
rz(-1.4830228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91325703) q[2];
sx q[2];
rz(-1.5568638) q[2];
sx q[2];
rz(3.081591) q[2];
rz(-1.2126728) q[3];
sx q[3];
rz(-2.4433177) q[3];
sx q[3];
rz(2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(2.5643964) q[0];
rz(-0.82950854) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(0.83782354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55623193) q[0];
sx q[0];
rz(-2.5149226) q[0];
sx q[0];
rz(0.13154948) q[0];
x q[1];
rz(-3.0810551) q[2];
sx q[2];
rz(-0.98145393) q[2];
sx q[2];
rz(-0.72414393) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41552222) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(0.46393053) q[1];
rz(0.44562037) q[3];
sx q[3];
rz(-0.94404781) q[3];
sx q[3];
rz(-2.8901951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(-1.2238097) q[2];
rz(-0.4661679) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849628) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(-3.0404941) q[0];
rz(2.7283607) q[1];
sx q[1];
rz(-0.44094545) q[1];
sx q[1];
rz(-2.0535927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4306385) q[0];
sx q[0];
rz(-2.0530409) q[0];
sx q[0];
rz(-2.7706465) q[0];
rz(-pi) q[1];
rz(-2.3067683) q[2];
sx q[2];
rz(-1.2918279) q[2];
sx q[2];
rz(2.5503412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9718379) q[1];
sx q[1];
rz(-2.7248451) q[1];
sx q[1];
rz(-0.93081148) q[1];
x q[2];
rz(2.4721778) q[3];
sx q[3];
rz(-1.6261887) q[3];
sx q[3];
rz(0.15633571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2130412) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(1.586033) q[2];
rz(-1.0726311) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17305408) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(-1.7065077) q[0];
rz(0.75812078) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(0.10163669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8827688) q[0];
sx q[0];
rz(-2.2166435) q[0];
sx q[0];
rz(2.1749138) q[0];
rz(-pi) q[1];
rz(0.91288699) q[2];
sx q[2];
rz(-0.55333558) q[2];
sx q[2];
rz(-1.1007512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58858777) q[1];
sx q[1];
rz(-1.2887452) q[1];
sx q[1];
rz(-1.4600919) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2472081) q[3];
sx q[3];
rz(-2.3830288) q[3];
sx q[3];
rz(2.0653084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3397843) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(1.3327117) q[2];
rz(2.0594275) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(-1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448755) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(2.3785059) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(0.9185763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42264274) q[0];
sx q[0];
rz(-3.115359) q[0];
sx q[0];
rz(1.8842741) q[0];
rz(2.1821676) q[2];
sx q[2];
rz(-2.583873) q[2];
sx q[2];
rz(-3.0447247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9715226) q[1];
sx q[1];
rz(-0.60107175) q[1];
sx q[1];
rz(-0.55723377) q[1];
x q[2];
rz(1.4537471) q[3];
sx q[3];
rz(-0.85005408) q[3];
sx q[3];
rz(0.085214867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(0.43583885) q[2];
rz(-0.11442746) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(-2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33335394) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(2.5575141) q[0];
rz(2.7059879) q[1];
sx q[1];
rz(-1.6807115) q[1];
sx q[1];
rz(1.6216507) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0348957) q[0];
sx q[0];
rz(-0.55904065) q[0];
sx q[0];
rz(-0.99026545) q[0];
rz(-pi) q[1];
rz(-2.9233515) q[2];
sx q[2];
rz(-1.1704966) q[2];
sx q[2];
rz(1.2117653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3518924) q[1];
sx q[1];
rz(-1.580211) q[1];
sx q[1];
rz(2.9355551) q[1];
x q[2];
rz(0.55840839) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(1.4608135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.074177563) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(2.0254859) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24935687) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(0.78999162) q[0];
rz(0.063591592) q[1];
sx q[1];
rz(-2.5345232) q[1];
sx q[1];
rz(-0.44073179) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249464) q[0];
sx q[0];
rz(-1.2727203) q[0];
sx q[0];
rz(-1.9008725) q[0];
rz(-2.47119) q[2];
sx q[2];
rz(-2.5385661) q[2];
sx q[2];
rz(-2.5489901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3778119) q[1];
sx q[1];
rz(-1.0965938) q[1];
sx q[1];
rz(-2.9567316) q[1];
rz(1.3578254) q[3];
sx q[3];
rz(-2.9479288) q[3];
sx q[3];
rz(1.97399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26238394) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(-0.66514307) q[2];
rz(2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36122286) q[0];
sx q[0];
rz(-2.4936115) q[0];
sx q[0];
rz(-1.004647) q[0];
rz(-1.7338344) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
sx q[1];
rz(1.2900603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72899517) q[0];
sx q[0];
rz(-0.83461232) q[0];
sx q[0];
rz(0.33552977) q[0];
rz(-pi) q[1];
rz(-1.8632221) q[2];
sx q[2];
rz(-2.5212581) q[2];
sx q[2];
rz(0.32493704) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5764783) q[1];
sx q[1];
rz(-1.8449515) q[1];
sx q[1];
rz(1.7903916) q[1];
rz(-pi) q[2];
rz(-1.3195924) q[3];
sx q[3];
rz(-2.525794) q[3];
sx q[3];
rz(-1.766966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0066234) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(2.8912365) q[2];
rz(2.1641459) q[3];
sx q[3];
rz(-1.2075295) q[3];
sx q[3];
rz(1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(0.74465887) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(-0.087564115) q[2];
sx q[2];
rz(-2.9976805) q[2];
sx q[2];
rz(1.9981071) q[2];
rz(-2.70784) q[3];
sx q[3];
rz(-1.1137782) q[3];
sx q[3];
rz(2.6657226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
