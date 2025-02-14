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
rz(-1.4327383) q[0];
sx q[0];
rz(-0.32810768) q[0];
sx q[0];
rz(-2.3018667) q[0];
rz(1.4108763) q[1];
sx q[1];
rz(-1.5411935) q[1];
sx q[1];
rz(0.76959258) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037207863) q[0];
sx q[0];
rz(-1.6055371) q[0];
sx q[0];
rz(-1.6503235) q[0];
rz(-0.02564557) q[2];
sx q[2];
rz(-0.8249976) q[2];
sx q[2];
rz(-2.1767669) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.283764) q[1];
sx q[1];
rz(-2.789838) q[1];
sx q[1];
rz(-1.5726388) q[1];
rz(-pi) q[2];
rz(-2.493164) q[3];
sx q[3];
rz(-2.3852013) q[3];
sx q[3];
rz(1.7738916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8771693) q[2];
sx q[2];
rz(-1.8921655) q[2];
sx q[2];
rz(1.1733615) q[2];
rz(2.7133283) q[3];
sx q[3];
rz(-2.5731125) q[3];
sx q[3];
rz(2.4885524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1332755) q[0];
sx q[0];
rz(-2.7011217) q[0];
sx q[0];
rz(-1.0622729) q[0];
rz(-1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(2.0468457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65829458) q[0];
sx q[0];
rz(-1.5196557) q[0];
sx q[0];
rz(-1.8258894) q[0];
rz(-pi) q[1];
rz(1.0350448) q[2];
sx q[2];
rz(-2.3887815) q[2];
sx q[2];
rz(0.24487803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1909263) q[1];
sx q[1];
rz(-2.6537173) q[1];
sx q[1];
rz(2.4937954) q[1];
rz(1.7844438) q[3];
sx q[3];
rz(-2.3348138) q[3];
sx q[3];
rz(0.90863228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4636479) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(2.2410683) q[2];
rz(-2.1053704) q[3];
sx q[3];
rz(-3.0039054) q[3];
sx q[3];
rz(-3.1305967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4765332) q[0];
sx q[0];
rz(-2.0158975) q[0];
sx q[0];
rz(1.4680468) q[0];
rz(-0.92848575) q[1];
sx q[1];
rz(-2.1378345) q[1];
sx q[1];
rz(-1.4220062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6674468) q[0];
sx q[0];
rz(-0.26413879) q[0];
sx q[0];
rz(2.788782) q[0];
rz(-pi) q[1];
rz(-0.81776889) q[2];
sx q[2];
rz(-2.1085133) q[2];
sx q[2];
rz(2.0995215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3513438) q[1];
sx q[1];
rz(-0.66748744) q[1];
sx q[1];
rz(0.60524596) q[1];
x q[2];
rz(1.8551781) q[3];
sx q[3];
rz(-1.1613628) q[3];
sx q[3];
rz(-2.6961498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15098393) q[2];
sx q[2];
rz(-1.432212) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(3.0153583) q[3];
sx q[3];
rz(-1.1773959) q[3];
sx q[3];
rz(1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6821297) q[0];
sx q[0];
rz(-1.4012902) q[0];
sx q[0];
rz(1.5392186) q[0];
rz(-1.0914717) q[1];
sx q[1];
rz(-0.83408728) q[1];
sx q[1];
rz(-1.9713255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6495958) q[0];
sx q[0];
rz(-1.2536712) q[0];
sx q[0];
rz(-1.1742422) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2762047) q[2];
sx q[2];
rz(-2.3995993) q[2];
sx q[2];
rz(-2.2592253) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2728426) q[1];
sx q[1];
rz(-0.12050546) q[1];
sx q[1];
rz(-2.4122448) q[1];
rz(2.6319632) q[3];
sx q[3];
rz(-0.42181236) q[3];
sx q[3];
rz(3.075656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22152659) q[2];
sx q[2];
rz(-0.92851323) q[2];
sx q[2];
rz(-3.0583734) q[2];
rz(2.4890238) q[3];
sx q[3];
rz(-3.1196399) q[3];
sx q[3];
rz(-0.86161247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49587747) q[0];
sx q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(-2.9448217) q[0];
rz(1.1605877) q[1];
sx q[1];
rz(-2.6996758) q[1];
sx q[1];
rz(-1.388185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3236448) q[0];
sx q[0];
rz(-2.9716688) q[0];
sx q[0];
rz(0.0090881149) q[0];
x q[1];
rz(0.20241995) q[2];
sx q[2];
rz(-1.7192063) q[2];
sx q[2];
rz(2.9925516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3922351) q[1];
sx q[1];
rz(-2.2171786) q[1];
sx q[1];
rz(2.9254167) q[1];
x q[2];
rz(1.9742244) q[3];
sx q[3];
rz(-2.3711172) q[3];
sx q[3];
rz(-1.4899173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56472003) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(1.8604856) q[2];
rz(0.34230226) q[3];
sx q[3];
rz(-3.0908995) q[3];
sx q[3];
rz(-2.2127693) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-1.0973278) q[0];
sx q[0];
rz(-1.0327562) q[0];
rz(0.67689854) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(-2.3597609) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1843518) q[0];
sx q[0];
rz(-0.68999422) q[0];
sx q[0];
rz(2.5646567) q[0];
rz(-2.135649) q[2];
sx q[2];
rz(-0.060533591) q[2];
sx q[2];
rz(-2.6356027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62877733) q[1];
sx q[1];
rz(-2.1330962) q[1];
sx q[1];
rz(-1.9356711) q[1];
rz(-2.8036814) q[3];
sx q[3];
rz(-1.5301782) q[3];
sx q[3];
rz(3.0062243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8402164) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(0.8737348) q[2];
rz(2.3524763) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(-2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18272045) q[0];
sx q[0];
rz(-0.046165753) q[0];
sx q[0];
rz(-0.090855457) q[0];
rz(0.60984045) q[1];
sx q[1];
rz(-1.5900541) q[1];
sx q[1];
rz(0.59744936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2029363) q[0];
sx q[0];
rz(-2.9466726) q[0];
sx q[0];
rz(-0.7233497) q[0];
x q[1];
rz(1.1979073) q[2];
sx q[2];
rz(-1.882236) q[2];
sx q[2];
rz(-0.58538891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3143464) q[1];
sx q[1];
rz(-1.1097849) q[1];
sx q[1];
rz(0.48996144) q[1];
rz(-1.1511943) q[3];
sx q[3];
rz(-1.6053995) q[3];
sx q[3];
rz(-2.7945003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7647543) q[2];
sx q[2];
rz(-0.54543269) q[2];
sx q[2];
rz(-3.0010014) q[2];
rz(-1.2815255) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-2.6144821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.918688) q[0];
sx q[0];
rz(-2.1253026) q[0];
sx q[0];
rz(-1.5284398) q[0];
rz(-2.3279922) q[1];
sx q[1];
rz(-0.10491144) q[1];
sx q[1];
rz(2.334107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19014004) q[0];
sx q[0];
rz(-1.7684816) q[0];
sx q[0];
rz(2.9485767) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1296223) q[2];
sx q[2];
rz(-1.1309012) q[2];
sx q[2];
rz(0.18378809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85881104) q[1];
sx q[1];
rz(-0.75560299) q[1];
sx q[1];
rz(-1.0048466) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2178879) q[3];
sx q[3];
rz(-2.4974303) q[3];
sx q[3];
rz(1.3180863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3537102) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(0.14459571) q[2];
rz(-2.0374129) q[3];
sx q[3];
rz(-2.927533) q[3];
sx q[3];
rz(-2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5398194) q[0];
sx q[0];
rz(-1.1310534) q[0];
sx q[0];
rz(0.55651504) q[0];
rz(-2.3717608) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(2.6803023) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249864) q[0];
sx q[0];
rz(-1.6184153) q[0];
sx q[0];
rz(-2.141593) q[0];
rz(-pi) q[1];
rz(0.67619063) q[2];
sx q[2];
rz(-0.95417038) q[2];
sx q[2];
rz(0.21756324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53992311) q[1];
sx q[1];
rz(-2.2846672) q[1];
sx q[1];
rz(1.9636092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57188028) q[3];
sx q[3];
rz(-1.9000305) q[3];
sx q[3];
rz(-2.5155932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7028659) q[2];
sx q[2];
rz(-2.8361969) q[2];
sx q[2];
rz(0.76496441) q[2];
rz(-1.2139828) q[3];
sx q[3];
rz(-2.5524804) q[3];
sx q[3];
rz(2.7629619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42534378) q[0];
sx q[0];
rz(-3.0916164) q[0];
sx q[0];
rz(-2.7783527) q[0];
rz(1.4092457) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(2.9149616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57412452) q[0];
sx q[0];
rz(-3.0063445) q[0];
sx q[0];
rz(0.02217205) q[0];
rz(0.93904597) q[2];
sx q[2];
rz(-0.9264731) q[2];
sx q[2];
rz(0.34560386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.781023) q[1];
sx q[1];
rz(-1.2464355) q[1];
sx q[1];
rz(0.85807577) q[1];
x q[2];
rz(-2.8037386) q[3];
sx q[3];
rz(-1.5607335) q[3];
sx q[3];
rz(-2.4110766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54214415) q[2];
sx q[2];
rz(-2.6207974) q[2];
sx q[2];
rz(0.64860541) q[2];
rz(3.0829698) q[3];
sx q[3];
rz(-0.62871814) q[3];
sx q[3];
rz(0.64529836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7764353) q[0];
sx q[0];
rz(-1.2189652) q[0];
sx q[0];
rz(-0.49268876) q[0];
rz(2.0877214) q[1];
sx q[1];
rz(-0.55640472) q[1];
sx q[1];
rz(0.41592204) q[1];
rz(-1.3345759) q[2];
sx q[2];
rz(-1.5470355) q[2];
sx q[2];
rz(2.7337337) q[2];
rz(-2.2661187) q[3];
sx q[3];
rz(-1.8322104) q[3];
sx q[3];
rz(-1.9598243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
