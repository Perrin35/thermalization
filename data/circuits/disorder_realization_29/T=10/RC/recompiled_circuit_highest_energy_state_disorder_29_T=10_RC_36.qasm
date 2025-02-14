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
rz(0.88275498) q[0];
sx q[0];
rz(-2.9967699) q[0];
sx q[0];
rz(0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061926024) q[0];
sx q[0];
rz(-0.57625895) q[0];
sx q[0];
rz(2.1012596) q[0];
rz(0.72291908) q[2];
sx q[2];
rz(-1.5462048) q[2];
sx q[2];
rz(1.493177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4262095) q[1];
sx q[1];
rz(-2.8236964) q[1];
sx q[1];
rz(0.23178394) q[1];
x q[2];
rz(3.0823067) q[3];
sx q[3];
rz(-0.93794696) q[3];
sx q[3];
rz(0.92801731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90193191) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(-2.530976) q[2];
rz(0.88879746) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438542) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(1.4854206) q[1];
sx q[1];
rz(-1.679436) q[1];
sx q[1];
rz(-0.86404538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9784911) q[0];
sx q[0];
rz(-1.3184883) q[0];
sx q[0];
rz(3.0358394) q[0];
x q[1];
rz(2.1435166) q[2];
sx q[2];
rz(-0.53600271) q[2];
sx q[2];
rz(2.0886476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5031918) q[1];
sx q[1];
rz(-1.4309037) q[1];
sx q[1];
rz(1.3963455) q[1];
x q[2];
rz(-0.95324272) q[3];
sx q[3];
rz(-0.88897486) q[3];
sx q[3];
rz(-0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(0.19134276) q[2];
rz(-2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3092344) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(2.7171296) q[0];
rz(-1.3407432) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(-2.4901966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021928) q[0];
sx q[0];
rz(-2.9777973) q[0];
sx q[0];
rz(1.6144362) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38369757) q[2];
sx q[2];
rz(-1.7406775) q[2];
sx q[2];
rz(0.023141247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60247682) q[1];
sx q[1];
rz(-1.0751564) q[1];
sx q[1];
rz(1.3631166) q[1];
x q[2];
rz(-1.797514) q[3];
sx q[3];
rz(-1.7245088) q[3];
sx q[3];
rz(-2.8938229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-2.885163) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7032787) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(2.1412204) q[0];
rz(1.4601624) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.9283074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98836658) q[0];
sx q[0];
rz(-1.8509764) q[0];
sx q[0];
rz(3.0154254) q[0];
rz(-pi) q[1];
x q[1];
rz(1.103066) q[2];
sx q[2];
rz(-2.2437895) q[2];
sx q[2];
rz(1.8813934) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4459232) q[1];
sx q[1];
rz(-1.4497821) q[1];
sx q[1];
rz(2.8112429) q[1];
rz(3.0509316) q[3];
sx q[3];
rz(-2.1075776) q[3];
sx q[3];
rz(1.2556374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(-0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475912) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(2.3498348) q[0];
rz(1.1119615) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70671088) q[0];
sx q[0];
rz(-1.1010873) q[0];
sx q[0];
rz(-0.89969866) q[0];
rz(0.15891517) q[2];
sx q[2];
rz(-2.5549485) q[2];
sx q[2];
rz(1.9288837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89692749) q[1];
sx q[1];
rz(-1.7411971) q[1];
sx q[1];
rz(-0.076431304) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8800635) q[3];
sx q[3];
rz(-0.59854186) q[3];
sx q[3];
rz(-0.33613294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2935334) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(2.000957) q[2];
rz(3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(-0.003224592) q[0];
rz(-0.047317304) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.3501732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4989717) q[0];
sx q[0];
rz(-1.823171) q[0];
sx q[0];
rz(-3.0809771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57729849) q[2];
sx q[2];
rz(-2.1566628) q[2];
sx q[2];
rz(-2.455204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9474831) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(-2.4653788) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7914823) q[3];
sx q[3];
rz(-1.6005362) q[3];
sx q[3];
rz(-1.1988175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99001592) q[2];
sx q[2];
rz(-2.9559957) q[2];
sx q[2];
rz(3.0604176) q[2];
rz(0.89933991) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.8544633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154295) q[0];
sx q[0];
rz(-1.8303215) q[0];
sx q[0];
rz(-2.6370866) q[0];
rz(0.99507487) q[1];
sx q[1];
rz(-2.1213396) q[1];
sx q[1];
rz(0.02034932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824469) q[0];
sx q[0];
rz(-1.2070281) q[0];
sx q[0];
rz(3.0412741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3624914) q[2];
sx q[2];
rz(-1.6951188) q[2];
sx q[2];
rz(-0.97036874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.523218) q[1];
sx q[1];
rz(-1.8357539) q[1];
sx q[1];
rz(-2.9441903) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0929111) q[3];
sx q[3];
rz(-1.4845856) q[3];
sx q[3];
rz(-1.3000559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67123479) q[2];
sx q[2];
rz(-1.3037668) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(1.5291519) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(-3.0488739) q[3];
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
rz(-2.2328211) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(2.0594647) q[0];
rz(-0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(-0.94295162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0531702) q[0];
sx q[0];
rz(-0.57500792) q[0];
sx q[0];
rz(-2.1413598) q[0];
rz(-2.2029938) q[2];
sx q[2];
rz(-1.7603121) q[2];
sx q[2];
rz(-1.0908443) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3228938) q[1];
sx q[1];
rz(-2.2320691) q[1];
sx q[1];
rz(-2.3803007) q[1];
x q[2];
rz(1.8887159) q[3];
sx q[3];
rz(-0.32264454) q[3];
sx q[3];
rz(1.2979729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97041398) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(1.9937493) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73087937) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-0.10243375) q[0];
rz(-0.61406413) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(-2.3238497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2376643) q[0];
sx q[0];
rz(-2.2127164) q[0];
sx q[0];
rz(-0.68591811) q[0];
rz(1.4253584) q[2];
sx q[2];
rz(-0.67611968) q[2];
sx q[2];
rz(-2.6349677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39108155) q[1];
sx q[1];
rz(-2.4121248) q[1];
sx q[1];
rz(1.5991648) q[1];
rz(2.2028186) q[3];
sx q[3];
rz(-1.9524487) q[3];
sx q[3];
rz(1.3248688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3488591) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(-0.46401986) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87026507) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(0.21367167) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(-0.22458354) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29447039) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(-0.76089528) q[0];
rz(-2.6170391) q[2];
sx q[2];
rz(-2.4894425) q[2];
sx q[2];
rz(0.99936501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.869369) q[1];
sx q[1];
rz(-1.0450796) q[1];
sx q[1];
rz(1.438526) q[1];
x q[2];
rz(2.0430327) q[3];
sx q[3];
rz(-2.3655112) q[3];
sx q[3];
rz(-0.14296338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97682041) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(-0.15178794) q[2];
rz(-3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0634154) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(2.0582485) q[1];
sx q[1];
rz(-1.5647519) q[1];
sx q[1];
rz(-1.5819989) q[1];
rz(0.2486817) q[2];
sx q[2];
rz(-2.9338825) q[2];
sx q[2];
rz(2.6354229) q[2];
rz(-0.32947551) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
