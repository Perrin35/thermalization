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
rz(-0.7769146) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033117753) q[0];
sx q[0];
rz(-0.43370789) q[0];
sx q[0];
rz(-1.2483622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5778779) q[2];
sx q[2];
rz(-1.6275121) q[2];
sx q[2];
rz(-0.67981718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0529671) q[1];
sx q[1];
rz(-1.8421116) q[1];
sx q[1];
rz(-1.1295098) q[1];
x q[2];
rz(0.99774811) q[3];
sx q[3];
rz(-1.5978509) q[3];
sx q[3];
rz(1.7535576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50600791) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(2.9721416) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(-1.0725526) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605543) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(3.0637528) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(0.30581623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77340375) q[0];
sx q[0];
rz(-1.0596674) q[0];
sx q[0];
rz(2.9298733) q[0];
rz(-pi) q[1];
rz(-0.14350899) q[2];
sx q[2];
rz(-0.75703758) q[2];
sx q[2];
rz(-0.17998634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9287323) q[1];
sx q[1];
rz(-1.1747735) q[1];
sx q[1];
rz(-0.81800445) q[1];
x q[2];
rz(-1.5068077) q[3];
sx q[3];
rz(-1.5100237) q[3];
sx q[3];
rz(1.2343386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1514312) q[2];
sx q[2];
rz(-1.7502681) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-3.0691872) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(-0.15637583) q[0];
rz(3.1001672) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(-1.5511537) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32937059) q[0];
sx q[0];
rz(-2.5539264) q[0];
sx q[0];
rz(2.7633322) q[0];
x q[1];
rz(-2.4855108) q[2];
sx q[2];
rz(-1.6245884) q[2];
sx q[2];
rz(-1.0840067) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2385905) q[1];
sx q[1];
rz(-2.6005473) q[1];
sx q[1];
rz(0.21958406) q[1];
rz(0.87167344) q[3];
sx q[3];
rz(-1.0314042) q[3];
sx q[3];
rz(1.2076857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40858832) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(-0.020922529) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9631831) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(0.43854976) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(-2.8964892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49798508) q[0];
sx q[0];
rz(-1.4873355) q[0];
sx q[0];
rz(2.2913234) q[0];
rz(-2.7736362) q[2];
sx q[2];
rz(-1.4237671) q[2];
sx q[2];
rz(-1.2556751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4192701) q[1];
sx q[1];
rz(-2.5825204) q[1];
sx q[1];
rz(-1.8354227) q[1];
rz(-pi) q[2];
rz(-2.80527) q[3];
sx q[3];
rz(-0.78236474) q[3];
sx q[3];
rz(-0.7366283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.601292) q[2];
sx q[2];
rz(-2.7056594) q[2];
sx q[2];
rz(-0.33176804) q[2];
rz(-0.48745421) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8496534) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(2.3680903) q[0];
rz(-1.9603112) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(-1.7519417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7252325) q[0];
sx q[0];
rz(-1.2073887) q[0];
sx q[0];
rz(-2.711722) q[0];
x q[1];
rz(0.59771363) q[2];
sx q[2];
rz(-2.1897912) q[2];
sx q[2];
rz(-0.60475443) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2323245) q[1];
sx q[1];
rz(-1.0977912) q[1];
sx q[1];
rz(-1.6632516) q[1];
rz(-pi) q[2];
rz(1.3895274) q[3];
sx q[3];
rz(-1.8049523) q[3];
sx q[3];
rz(-0.93269809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9224077) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(2.3310272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(-0.26350185) q[0];
rz(1.1031411) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(-0.37364328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140681) q[0];
sx q[0];
rz(-1.6434945) q[0];
sx q[0];
rz(3.0811653) q[0];
x q[1];
rz(0.12219001) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(3.0360589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9571324) q[1];
sx q[1];
rz(-1.2533979) q[1];
sx q[1];
rz(-0.98021345) q[1];
x q[2];
rz(2.13426) q[3];
sx q[3];
rz(-2.084199) q[3];
sx q[3];
rz(2.8443984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0282447) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(2.587758) q[2];
rz(-1.7438186) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(-2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(-0.30034932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25324437) q[0];
sx q[0];
rz(-1.7333366) q[0];
sx q[0];
rz(-1.4534611) q[0];
x q[1];
rz(1.0763758) q[2];
sx q[2];
rz(-2.2110614) q[2];
sx q[2];
rz(-0.51039052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19076599) q[1];
sx q[1];
rz(-0.015282282) q[1];
sx q[1];
rz(-1.5453668) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85791479) q[3];
sx q[3];
rz(-1.5893717) q[3];
sx q[3];
rz(-1.9481079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.131669) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(-0.24492502) q[2];
rz(-0.51982546) q[3];
sx q[3];
rz(-2.2782785) q[3];
sx q[3];
rz(-2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7898665) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(1.1630195) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(-1.012872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236008) q[0];
sx q[0];
rz(-2.3396684) q[0];
sx q[0];
rz(0.84419925) q[0];
rz(-pi) q[1];
rz(-1.0032907) q[2];
sx q[2];
rz(-1.419628) q[2];
sx q[2];
rz(1.0790107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.098274577) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(1.2374452) q[1];
rz(-1.1539641) q[3];
sx q[3];
rz(-2.2139021) q[3];
sx q[3];
rz(-2.7288306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11671994) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087273) q[0];
sx q[0];
rz(-2.9948586) q[0];
sx q[0];
rz(3.1291381) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(0.27997231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075172193) q[0];
sx q[0];
rz(-3.018258) q[0];
sx q[0];
rz(-1.9955817) q[0];
rz(-0.32692636) q[2];
sx q[2];
rz(-1.1571615) q[2];
sx q[2];
rz(1.3437611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1713531) q[1];
sx q[1];
rz(-2.2917124) q[1];
sx q[1];
rz(1.156714) q[1];
rz(0.98484184) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(-2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81194699) q[2];
sx q[2];
rz(-2.3880366) q[2];
sx q[2];
rz(-0.21198708) q[2];
rz(-2.3181465) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(0.25920355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9987746) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(0.69277358) q[0];
rz(0.57299262) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(0.43100345) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26769022) q[0];
sx q[0];
rz(-0.55492102) q[0];
sx q[0];
rz(0.11124994) q[0];
rz(1.1262604) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(-2.0972507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17934701) q[1];
sx q[1];
rz(-1.822346) q[1];
sx q[1];
rz(-2.4169282) q[1];
rz(-pi) q[2];
rz(-2.1206843) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(-2.6481215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(-0.46960056) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(-0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(2.8259059) q[2];
sx q[2];
rz(-1.9236947) q[2];
sx q[2];
rz(-2.535939) q[2];
rz(-2.2297494) q[3];
sx q[3];
rz(-1.5921436) q[3];
sx q[3];
rz(2.5719503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
