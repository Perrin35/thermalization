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
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(2.3646781) q[0];
rz(-1.7493526) q[1];
sx q[1];
rz(-1.8267781) q[1];
sx q[1];
rz(-2.1652752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1084749) q[0];
sx q[0];
rz(-2.7078848) q[0];
sx q[0];
rz(-1.2483622) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10586057) q[2];
sx q[2];
rz(-2.5753394) q[2];
sx q[2];
rz(-0.98041269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9982517) q[1];
sx q[1];
rz(-0.51330459) q[1];
sx q[1];
rz(-2.1480888) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99774811) q[3];
sx q[3];
rz(-1.5978509) q[3];
sx q[3];
rz(-1.7535576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50600791) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(0.40679833) q[2];
rz(-2.9721416) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(-2.0690401) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88103831) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(0.30581623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394442) q[0];
sx q[0];
rz(-1.3864707) q[0];
sx q[0];
rz(2.091616) q[0];
x q[1];
rz(-0.14350899) q[2];
sx q[2];
rz(-2.3845551) q[2];
sx q[2];
rz(-2.9616063) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0111573) q[1];
sx q[1];
rz(-2.2533335) q[1];
sx q[1];
rz(-2.621317) q[1];
rz(-pi) q[2];
rz(2.3313794) q[3];
sx q[3];
rz(-3.0533724) q[3];
sx q[3];
rz(2.719413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1514312) q[2];
sx q[2];
rz(-1.7502681) q[2];
sx q[2];
rz(2.4988417) q[2];
rz(3.0691872) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0533503) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(0.15637583) q[0];
rz(3.1001672) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(-1.590439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0254733) q[0];
sx q[0];
rz(-2.1120694) q[0];
sx q[0];
rz(1.812029) q[0];
rz(-pi) q[1];
rz(-0.65608187) q[2];
sx q[2];
rz(-1.5170043) q[2];
sx q[2];
rz(2.057586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.620365) q[1];
sx q[1];
rz(-1.4583734) q[1];
sx q[1];
rz(2.6112154) q[1];
rz(-0.87167344) q[3];
sx q[3];
rz(-2.1101885) q[3];
sx q[3];
rz(1.2076857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7330043) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(-3.1206701) q[2];
rz(-0.18713348) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9631831) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(-0.43854976) q[0];
rz(1.5248388) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(2.8964892) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995718) q[0];
sx q[0];
rz(-0.85332131) q[0];
sx q[0];
rz(3.0307253) q[0];
rz(2.7510277) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(0.67827889) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72232258) q[1];
sx q[1];
rz(-2.5825204) q[1];
sx q[1];
rz(1.8354227) q[1];
rz(-pi) q[2];
rz(0.33632261) q[3];
sx q[3];
rz(-2.3592279) q[3];
sx q[3];
rz(-2.4049644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54030067) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(2.8098246) q[2];
rz(-2.6541384) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8496534) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(2.3680903) q[0];
rz(1.1812814) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(-1.7519417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4163602) q[0];
sx q[0];
rz(-1.9342039) q[0];
sx q[0];
rz(-0.4298707) q[0];
rz(-pi) q[1];
rz(0.90221407) q[2];
sx q[2];
rz(-0.83219516) q[2];
sx q[2];
rz(2.8813643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4378499) q[1];
sx q[1];
rz(-1.653076) q[1];
sx q[1];
rz(-0.47474307) q[1];
rz(0.23791194) q[3];
sx q[3];
rz(-1.7470659) q[3];
sx q[3];
rz(-2.545994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(0.47214559) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.8483714) q[3];
sx q[3];
rz(-2.3310272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(-1.1031411) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(2.7679494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81955302) q[0];
sx q[0];
rz(-3.0470938) q[0];
sx q[0];
rz(-2.2631133) q[0];
rz(-pi) q[1];
rz(3.0194026) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(-3.0360589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3195575) q[1];
sx q[1];
rz(-2.4801835) q[1];
sx q[1];
rz(1.0378077) q[1];
x q[2];
rz(-2.3832537) q[3];
sx q[3];
rz(-2.3985574) q[3];
sx q[3];
rz(-0.61279994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11334795) q[2];
sx q[2];
rz(-0.17013203) q[2];
sx q[2];
rz(0.55383468) q[2];
rz(1.7438186) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5572307) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.2297909) q[0];
rz(-2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(0.30034932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984788) q[0];
sx q[0];
rz(-1.686578) q[0];
sx q[0];
rz(-0.16364574) q[0];
x q[1];
rz(-0.70234583) q[2];
sx q[2];
rz(-1.9611729) q[2];
sx q[2];
rz(2.3927488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19076599) q[1];
sx q[1];
rz(-0.015282282) q[1];
sx q[1];
rz(-1.5453668) q[1];
rz(-pi) q[2];
rz(-2.2836779) q[3];
sx q[3];
rz(-1.5893717) q[3];
sx q[3];
rz(-1.1934848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0099237) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(2.8966676) q[2];
rz(-2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(-2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7898665) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(1.9785731) q[0];
rz(-3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(1.012872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.430535) q[0];
sx q[0];
rz(-1.0037046) q[0];
sx q[0];
rz(2.5398769) q[0];
x q[1];
rz(1.8469454) q[2];
sx q[2];
rz(-2.5564402) q[2];
sx q[2];
rz(0.25979751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7808395) q[1];
sx q[1];
rz(-1.8855125) q[1];
sx q[1];
rz(0.34784045) q[1];
x q[2];
rz(-1.9876285) q[3];
sx q[3];
rz(-0.92769054) q[3];
sx q[3];
rz(-2.7288306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0248727) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(0.60574496) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(-0.34887031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054319687) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(3.1291381) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-0.27997231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6387647) q[0];
sx q[0];
rz(-1.6831213) q[0];
sx q[0];
rz(3.0905484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8146663) q[2];
sx q[2];
rz(-1.1571615) q[2];
sx q[2];
rz(-1.7978316) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1713531) q[1];
sx q[1];
rz(-2.2917124) q[1];
sx q[1];
rz(-1.9848787) q[1];
x q[2];
rz(1.9334698) q[3];
sx q[3];
rz(-1.8020523) q[3];
sx q[3];
rz(-1.4303007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3296457) q[2];
sx q[2];
rz(-2.3880366) q[2];
sx q[2];
rz(2.9296056) q[2];
rz(2.3181465) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(-2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.9987746) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(0.69277358) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(0.43100345) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-0.55492102) q[0];
sx q[0];
rz(-3.0303427) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1262604) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(1.044342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5332408) q[1];
sx q[1];
rz(-2.2679288) q[1];
sx q[1];
rz(-1.240154) q[1];
rz(-2.1206843) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(0.49347116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(2.4277021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8568759) q[0];
sx q[0];
rz(-1.4061883) q[0];
sx q[0];
rz(-1.0549369) q[0];
rz(-0.84125413) q[1];
sx q[1];
rz(-2.0396736) q[1];
sx q[1];
rz(3.094818) q[1];
rz(1.2011436) q[2];
sx q[2];
rz(-1.2751725) q[2];
sx q[2];
rz(-1.0775492) q[2];
rz(1.6056521) q[3];
sx q[3];
rz(-0.65924725) q[3];
sx q[3];
rz(1.02871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
