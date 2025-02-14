OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(-1.9174175) q[0];
sx q[0];
rz(1.8608215) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(0.26689902) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163942) q[0];
sx q[0];
rz(-0.13422899) q[0];
sx q[0];
rz(0.79948808) q[0];
x q[1];
rz(-2.884567) q[2];
sx q[2];
rz(-0.89867175) q[2];
sx q[2];
rz(2.4804403) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7952629) q[1];
sx q[1];
rz(-0.85546934) q[1];
sx q[1];
rz(2.1500509) q[1];
rz(-pi) q[2];
rz(0.21453114) q[3];
sx q[3];
rz(-0.45154587) q[3];
sx q[3];
rz(-0.2487693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12307564) q[2];
sx q[2];
rz(-2.1776431) q[2];
sx q[2];
rz(0.71259552) q[2];
rz(2.9739042) q[3];
sx q[3];
rz(-0.84961397) q[3];
sx q[3];
rz(0.4445506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0406822) q[0];
sx q[0];
rz(-1.7829144) q[0];
sx q[0];
rz(-1.6810625) q[0];
rz(1.4617317) q[1];
sx q[1];
rz(-1.0667543) q[1];
sx q[1];
rz(2.0251958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76976639) q[0];
sx q[0];
rz(-0.049749181) q[0];
sx q[0];
rz(-2.0905963) q[0];
rz(-pi) q[1];
rz(1.0149392) q[2];
sx q[2];
rz(-1.6995322) q[2];
sx q[2];
rz(-1.577108) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18552407) q[1];
sx q[1];
rz(-1.9466234) q[1];
sx q[1];
rz(1.8517428) q[1];
rz(-pi) q[2];
rz(-2.8380409) q[3];
sx q[3];
rz(-2.7139691) q[3];
sx q[3];
rz(2.7855765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98550335) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(-2.8893341) q[2];
rz(-1.0559399) q[3];
sx q[3];
rz(-0.85774937) q[3];
sx q[3];
rz(-1.0004388) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58241874) q[0];
sx q[0];
rz(-2.2370339) q[0];
sx q[0];
rz(2.3140123) q[0];
rz(-0.52455348) q[1];
sx q[1];
rz(-1.2453715) q[1];
sx q[1];
rz(0.96447271) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8480769) q[0];
sx q[0];
rz(-1.6381253) q[0];
sx q[0];
rz(0.99253722) q[0];
rz(-pi) q[1];
rz(-3.0542637) q[2];
sx q[2];
rz(-1.1967518) q[2];
sx q[2];
rz(2.5279074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71389329) q[1];
sx q[1];
rz(-1.5212219) q[1];
sx q[1];
rz(-2.5109641) q[1];
rz(-pi) q[2];
rz(-3.019624) q[3];
sx q[3];
rz(-0.44443529) q[3];
sx q[3];
rz(1.1201064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7437637) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(1.942983) q[2];
rz(-1.4926636) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(-0.57884136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32597932) q[0];
sx q[0];
rz(-2.9878243) q[0];
sx q[0];
rz(-1.9258668) q[0];
rz(-2.1513596) q[1];
sx q[1];
rz(-2.4001887) q[1];
sx q[1];
rz(0.56891099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7051681) q[0];
sx q[0];
rz(-1.2878083) q[0];
sx q[0];
rz(0.97375906) q[0];
rz(-pi) q[1];
rz(0.26428916) q[2];
sx q[2];
rz(-1.2131872) q[2];
sx q[2];
rz(2.8823095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4601656) q[1];
sx q[1];
rz(-2.6003464) q[1];
sx q[1];
rz(2.1591004) q[1];
rz(-pi) q[2];
rz(-0.74796933) q[3];
sx q[3];
rz(-2.0828205) q[3];
sx q[3];
rz(-1.0486163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8625921) q[2];
sx q[2];
rz(-2.180884) q[2];
sx q[2];
rz(2.847239) q[2];
rz(-2.477008) q[3];
sx q[3];
rz(-2.3817101) q[3];
sx q[3];
rz(1.6871281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6890474) q[0];
sx q[0];
rz(-1.713546) q[0];
sx q[0];
rz(2.8277165) q[0];
rz(2.2398056) q[1];
sx q[1];
rz(-1.7867463) q[1];
sx q[1];
rz(1.6759759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4340227) q[0];
sx q[0];
rz(-1.5965874) q[0];
sx q[0];
rz(-1.4517734) q[0];
x q[1];
rz(1.1526665) q[2];
sx q[2];
rz(-0.2731495) q[2];
sx q[2];
rz(-1.7371617) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82354083) q[1];
sx q[1];
rz(-2.7465804) q[1];
sx q[1];
rz(-0.1956296) q[1];
x q[2];
rz(3.0907822) q[3];
sx q[3];
rz(-0.63119315) q[3];
sx q[3];
rz(0.3785336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19447154) q[2];
sx q[2];
rz(-2.3919545) q[2];
sx q[2];
rz(2.0571902) q[2];
rz(-0.32367745) q[3];
sx q[3];
rz(-2.4249228) q[3];
sx q[3];
rz(0.22423854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913469) q[0];
sx q[0];
rz(-0.24557376) q[0];
sx q[0];
rz(-0.41159758) q[0];
rz(0.72548524) q[1];
sx q[1];
rz(-1.8975703) q[1];
sx q[1];
rz(-1.4010319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4832725) q[0];
sx q[0];
rz(-1.3764769) q[0];
sx q[0];
rz(-0.65388314) q[0];
rz(-pi) q[1];
rz(-2.4011924) q[2];
sx q[2];
rz(-1.4776738) q[2];
sx q[2];
rz(-1.4880866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7159375) q[1];
sx q[1];
rz(-1.2285985) q[1];
sx q[1];
rz(-1.2340496) q[1];
x q[2];
rz(-1.6468923) q[3];
sx q[3];
rz(-2.020936) q[3];
sx q[3];
rz(-1.744818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97658849) q[2];
sx q[2];
rz(-2.5133331) q[2];
sx q[2];
rz(-1.7151625) q[2];
rz(-0.12380883) q[3];
sx q[3];
rz(-1.0694458) q[3];
sx q[3];
rz(-2.6482705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8641149) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(1.3364828) q[0];
rz(-0.9388963) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(-0.28405651) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1079252) q[0];
sx q[0];
rz(-1.9104275) q[0];
sx q[0];
rz(2.4976298) q[0];
rz(-pi) q[1];
rz(-1.797126) q[2];
sx q[2];
rz(-1.7301736) q[2];
sx q[2];
rz(1.2633737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56834451) q[1];
sx q[1];
rz(-1.580984) q[1];
sx q[1];
rz(-1.6022026) q[1];
x q[2];
rz(1.7098996) q[3];
sx q[3];
rz(-1.3918607) q[3];
sx q[3];
rz(3.1070955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9629918) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(2.2006688) q[2];
rz(-3.1407147) q[3];
sx q[3];
rz(-1.2255171) q[3];
sx q[3];
rz(-0.0860478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0963928) q[0];
sx q[0];
rz(-0.89458507) q[0];
sx q[0];
rz(-2.9659502) q[0];
rz(1.2200217) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(2.2697935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052446121) q[0];
sx q[0];
rz(-1.8348357) q[0];
sx q[0];
rz(-2.9857078) q[0];
rz(0.7824509) q[2];
sx q[2];
rz(-1.8385213) q[2];
sx q[2];
rz(-0.29111448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3094191) q[1];
sx q[1];
rz(-1.4117575) q[1];
sx q[1];
rz(-2.5320903) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55671911) q[3];
sx q[3];
rz(-1.5079323) q[3];
sx q[3];
rz(-1.6281466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4153727) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(-2.1190763) q[2];
rz(2.7546049) q[3];
sx q[3];
rz(-1.756668) q[3];
sx q[3];
rz(-2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047394) q[0];
sx q[0];
rz(-1.1664718) q[0];
sx q[0];
rz(0.26300305) q[0];
rz(-2.8291342) q[1];
sx q[1];
rz(-1.0954906) q[1];
sx q[1];
rz(-1.1824664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1240886) q[0];
sx q[0];
rz(-2.5944983) q[0];
sx q[0];
rz(-1.9518243) q[0];
x q[1];
rz(0.55652852) q[2];
sx q[2];
rz(-0.91224837) q[2];
sx q[2];
rz(1.5040894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8856462) q[1];
sx q[1];
rz(-1.3722536) q[1];
sx q[1];
rz(0.37862733) q[1];
x q[2];
rz(0.83689383) q[3];
sx q[3];
rz(-0.7115127) q[3];
sx q[3];
rz(1.6435028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.075228127) q[2];
sx q[2];
rz(-1.3508537) q[2];
sx q[2];
rz(-0.24825516) q[2];
rz(-2.9634641) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(0.78290141) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42661509) q[0];
sx q[0];
rz(-0.39418945) q[0];
sx q[0];
rz(-2.07975) q[0];
rz(1.3007523) q[1];
sx q[1];
rz(-1.6164833) q[1];
sx q[1];
rz(-1.1311857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12515629) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(-0.50986503) q[0];
x q[1];
rz(1.1243368) q[2];
sx q[2];
rz(-2.0314337) q[2];
sx q[2];
rz(1.8370475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0319236) q[1];
sx q[1];
rz(-1.1732035) q[1];
sx q[1];
rz(-1.6709187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33343306) q[3];
sx q[3];
rz(-2.5134183) q[3];
sx q[3];
rz(-1.4905358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0608369) q[2];
sx q[2];
rz(-0.26641014) q[2];
sx q[2];
rz(2.5401435) q[2];
rz(1.0414177) q[3];
sx q[3];
rz(-1.6626549) q[3];
sx q[3];
rz(-1.3781594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.9926485) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(3.0638937) q[1];
sx q[1];
rz(-2.6014889) q[1];
sx q[1];
rz(-1.4059975) q[1];
rz(-0.97066047) q[2];
sx q[2];
rz(-1.9117336) q[2];
sx q[2];
rz(-2.5993549) q[2];
rz(-0.79176767) q[3];
sx q[3];
rz(-1.0825915) q[3];
sx q[3];
rz(-1.3655567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
