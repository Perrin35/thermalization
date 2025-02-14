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
rz(0.32938862) q[0];
sx q[0];
rz(-2.5100799) q[0];
sx q[0];
rz(0.00087498571) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(2.7963918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314197) q[0];
sx q[0];
rz(-1.6616482) q[0];
sx q[0];
rz(-3.1167555) q[0];
rz(-pi) q[1];
rz(2.7852374) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(-0.25549437) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4507323) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(2.2295206) q[1];
rz(1.6451938) q[3];
sx q[3];
rz(-1.5545903) q[3];
sx q[3];
rz(-2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-2.148707) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(-2.3670926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82371432) q[0];
sx q[0];
rz(-1.5597938) q[0];
sx q[0];
rz(-1.3925793) q[0];
rz(-pi) q[1];
rz(2.3341137) q[2];
sx q[2];
rz(-1.1801071) q[2];
sx q[2];
rz(1.1451858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5238873) q[1];
sx q[1];
rz(-1.3464217) q[1];
sx q[1];
rz(0.38274204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9162972) q[3];
sx q[3];
rz(-2.0301798) q[3];
sx q[3];
rz(0.53250203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.896686) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(3.0027622) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15181825) q[0];
sx q[0];
rz(-2.646137) q[0];
sx q[0];
rz(2.302343) q[0];
rz(-pi) q[1];
rz(-2.6723318) q[2];
sx q[2];
rz(-2.3692694) q[2];
sx q[2];
rz(-1.1132211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.638354) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(2.9364763) q[1];
x q[2];
rz(-2.5957727) q[3];
sx q[3];
rz(-2.8132317) q[3];
sx q[3];
rz(0.92121802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46131721) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(2.7834564) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0963652) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(-0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-0.92794424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5281677) q[0];
sx q[0];
rz(-1.6990464) q[0];
sx q[0];
rz(1.7665461) q[0];
rz(-pi) q[1];
rz(-2.5408542) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(-1.635765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18670652) q[1];
sx q[1];
rz(-2.2518603) q[1];
sx q[1];
rz(-1.7658556) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0870861) q[3];
sx q[3];
rz(-1.0060423) q[3];
sx q[3];
rz(1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(-2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(-1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(1.060932) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.2394989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0415619) q[0];
sx q[0];
rz(-1.8503555) q[0];
sx q[0];
rz(-0.97645219) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4475736) q[2];
sx q[2];
rz(-0.87022129) q[2];
sx q[2];
rz(-1.4365591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8722788) q[1];
sx q[1];
rz(-1.8511103) q[1];
sx q[1];
rz(2.6107236) q[1];
rz(-pi) q[2];
rz(-1.0723713) q[3];
sx q[3];
rz(-1.9349422) q[3];
sx q[3];
rz(-2.2143603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(-1.3266374) q[2];
rz(0.2462247) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-2.2897913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(-2.7958561) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-0.87127042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7428674) q[0];
sx q[0];
rz(-0.79116066) q[0];
sx q[0];
rz(-1.4732811) q[0];
rz(-0.75052336) q[2];
sx q[2];
rz(-1.125968) q[2];
sx q[2];
rz(-1.9409279) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6253736) q[1];
sx q[1];
rz(-1.1828848) q[1];
sx q[1];
rz(-0.078887786) q[1];
rz(1.096181) q[3];
sx q[3];
rz(-2.1283009) q[3];
sx q[3];
rz(-1.6826009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8295916) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(-1.568659) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0503814) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0480014) q[0];
sx q[0];
rz(-2.8107852) q[0];
sx q[0];
rz(0.76865102) q[0];
rz(-2.4270646) q[2];
sx q[2];
rz(-2.0863669) q[2];
sx q[2];
rz(2.553568) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70434785) q[1];
sx q[1];
rz(-1.7926551) q[1];
sx q[1];
rz(1.1997486) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0934392) q[3];
sx q[3];
rz(-0.94814903) q[3];
sx q[3];
rz(-2.4751055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1130134) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(2.8431559) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(2.723208) q[0];
rz(-2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(3.0072838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1100562) q[0];
sx q[0];
rz(-0.93288619) q[0];
sx q[0];
rz(2.4619034) q[0];
x q[1];
rz(1.1654794) q[2];
sx q[2];
rz(-1.3994872) q[2];
sx q[2];
rz(1.5660945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.011439104) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(0.18421872) q[1];
rz(-pi) q[2];
rz(2.8577096) q[3];
sx q[3];
rz(-2.0238658) q[3];
sx q[3];
rz(0.69597746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(1.0632473) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(-1.0412019) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5409656) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-1.0293915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5350069) q[0];
sx q[0];
rz(-0.86678737) q[0];
sx q[0];
rz(-0.042198472) q[0];
x q[1];
rz(1.5310982) q[2];
sx q[2];
rz(-1.2671114) q[2];
sx q[2];
rz(0.025321753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16949108) q[1];
sx q[1];
rz(-1.6147366) q[1];
sx q[1];
rz(1.340534) q[1];
x q[2];
rz(-1.6624032) q[3];
sx q[3];
rz(-1.1044958) q[3];
sx q[3];
rz(-0.3084076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(-0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(2.9990101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.044238) q[0];
sx q[0];
rz(-1.729106) q[0];
sx q[0];
rz(1.5990785) q[0];
x q[1];
rz(-2.768553) q[2];
sx q[2];
rz(-2.4786886) q[2];
sx q[2];
rz(0.66696862) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8293162) q[1];
sx q[1];
rz(-1.3061151) q[1];
sx q[1];
rz(0.97977248) q[1];
rz(-pi) q[2];
rz(2.1964588) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(1.6459203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(-1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0161229) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(0.62025537) q[2];
sx q[2];
rz(-0.56369416) q[2];
sx q[2];
rz(-0.77947215) q[2];
rz(-3.075243) q[3];
sx q[3];
rz(-1.2766311) q[3];
sx q[3];
rz(1.0286812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
