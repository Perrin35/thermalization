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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(-0.00087498571) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314197) q[0];
sx q[0];
rz(-1.6616482) q[0];
sx q[0];
rz(0.024837107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7183185) q[2];
sx q[2];
rz(-1.194724) q[2];
sx q[2];
rz(-3.0126115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7789014) q[1];
sx q[1];
rz(-0.76299113) q[1];
sx q[1];
rz(-0.94339006) q[1];
rz(1.3561068) q[3];
sx q[3];
rz(-3.0654538) q[3];
sx q[3];
rz(2.3977181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(-1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107373) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(2.3670926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607722) q[0];
sx q[0];
rz(-2.9630399) q[0];
sx q[0];
rz(1.6327842) q[0];
rz(0.51807816) q[2];
sx q[2];
rz(-0.87730125) q[2];
sx q[2];
rz(-2.3665646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0424287) q[1];
sx q[1];
rz(-1.9434669) q[1];
sx q[1];
rz(-1.3295688) q[1];
x q[2];
rz(1.9949739) q[3];
sx q[3];
rz(-0.50809233) q[3];
sx q[3];
rz(1.0095694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.896686) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(1.0150821) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-0.15750289) q[0];
rz(-2.1306254) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(3.0027622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75051266) q[0];
sx q[0];
rz(-1.8939928) q[0];
sx q[0];
rz(1.1884407) q[0];
rz(-0.71535297) q[2];
sx q[2];
rz(-1.2497447) q[2];
sx q[2];
rz(-2.3356444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0142747) q[1];
sx q[1];
rz(-1.7671314) q[1];
sx q[1];
rz(1.2731958) q[1];
rz(-pi) q[2];
rz(-0.54581996) q[3];
sx q[3];
rz(-0.32836093) q[3];
sx q[3];
rz(-2.2203746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(-2.7834564) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(1.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(0.3048234) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-2.2136484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6154895) q[0];
sx q[0];
rz(-2.9080221) q[0];
sx q[0];
rz(-2.156267) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60073845) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(1.5058277) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18670652) q[1];
sx q[1];
rz(-0.88973239) q[1];
sx q[1];
rz(1.7658556) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66189142) q[3];
sx q[3];
rz(-0.74569476) q[3];
sx q[3];
rz(0.5058561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(0.18816571) q[2];
rz(1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(-1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6998049) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(-1.2394989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0584377) q[0];
sx q[0];
rz(-0.64955901) q[0];
sx q[0];
rz(-2.044528) q[0];
rz(-pi) q[1];
rz(-2.4475736) q[2];
sx q[2];
rz(-0.87022129) q[2];
sx q[2];
rz(1.7050336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2693138) q[1];
sx q[1];
rz(-1.8511103) q[1];
sx q[1];
rz(-0.53086908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0692213) q[3];
sx q[3];
rz(-1.9349422) q[3];
sx q[3];
rz(2.2143603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061848) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-0.87127042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987253) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(1.6683116) q[0];
rz(-0.99314697) q[2];
sx q[2];
rz(-0.9075853) q[2];
sx q[2];
rz(0.011485966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(1.9597998) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61136758) q[3];
sx q[3];
rz(-1.9690367) q[3];
sx q[3];
rz(-0.37722019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8295916) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(0.87257067) q[2];
rz(-1.568659) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(-0.76164371) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0480014) q[0];
sx q[0];
rz(-2.8107852) q[0];
sx q[0];
rz(-0.76865102) q[0];
rz(-pi) q[1];
rz(-0.71304597) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(2.6756659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35173479) q[1];
sx q[1];
rz(-2.7119293) q[1];
sx q[1];
rz(-1.0142782) q[1];
x q[2];
rz(-1.6377444) q[3];
sx q[3];
rz(-2.5173325) q[3];
sx q[3];
rz(0.74893307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1130134) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(2.8431559) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(-2.723208) q[0];
rz(2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1100562) q[0];
sx q[0];
rz(-2.2087065) q[0];
sx q[0];
rz(-2.4619034) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1573345) q[2];
sx q[2];
rz(-0.43817876) q[2];
sx q[2];
rz(-0.38288051) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5843552) q[1];
sx q[1];
rz(-1.3882625) q[1];
sx q[1];
rz(1.4339158) q[1];
rz(1.1014492) q[3];
sx q[3];
rz(-1.3162287) q[3];
sx q[3];
rz(2.1397487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7316651) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(2.4995787) q[2];
rz(-2.0783453) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(-1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5409656) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(3.0294321) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(1.0293915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0784796) q[0];
sx q[0];
rz(-1.6029583) q[0];
sx q[0];
rz(-2.2752447) q[0];
rz(-pi) q[1];
rz(0.30390988) q[2];
sx q[2];
rz(-1.5329156) q[2];
sx q[2];
rz(-1.5335976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7299906) q[1];
sx q[1];
rz(-1.8008324) q[1];
sx q[1];
rz(-3.0964628) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9618046) q[3];
sx q[3];
rz(-0.4745634) q[3];
sx q[3];
rz(-3.0347412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42334291) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(-2.5293317) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(-1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(-2.9990101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6194942) q[0];
sx q[0];
rz(-1.5987248) q[0];
sx q[0];
rz(0.15837196) q[0];
rz(-2.768553) q[2];
sx q[2];
rz(-2.4786886) q[2];
sx q[2];
rz(0.66696862) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0568472) q[1];
sx q[1];
rz(-2.1386302) q[1];
sx q[1];
rz(-0.3155057) q[1];
x q[2];
rz(-2.1964588) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(1.4956724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3451781) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(-1.4043407) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1254697) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(1.2186981) q[2];
sx q[2];
rz(-1.1209956) q[2];
sx q[2];
rz(1.6605177) q[2];
rz(-0.066349647) q[3];
sx q[3];
rz(-1.8649615) q[3];
sx q[3];
rz(-2.1129114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
