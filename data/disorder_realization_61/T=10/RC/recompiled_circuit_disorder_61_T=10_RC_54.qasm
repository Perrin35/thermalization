OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(1.1336741) q[0];
sx q[0];
rz(7.8757474) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(2.5974098) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2092752) q[0];
sx q[0];
rz(-1.5341175) q[0];
sx q[0];
rz(3.0741312) q[0];
rz(-0.10138301) q[2];
sx q[2];
rz(-1.6147436) q[2];
sx q[2];
rz(-1.6909388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4197598) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(-2.1268197) q[1];
rz(-pi) q[2];
rz(2.5892026) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(-2.606483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(1.7791746) q[2];
rz(0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966999) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(1.171296) q[0];
rz(-0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(1.404095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(1.0147592) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7369466) q[2];
sx q[2];
rz(-1.0937905) q[2];
sx q[2];
rz(-2.0376861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5215056) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(2.8584245) q[1];
rz(-pi) q[2];
rz(-2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8319548) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(2.3220093) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-2.5779285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6113341) q[0];
sx q[0];
rz(-1.5000492) q[0];
sx q[0];
rz(2.5669839) q[0];
rz(-0.59994772) q[2];
sx q[2];
rz(-1.3834582) q[2];
sx q[2];
rz(-1.7768605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3445661) q[1];
sx q[1];
rz(-1.016942) q[1];
sx q[1];
rz(0.38573854) q[1];
rz(1.1590957) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(1.2702277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(-2.0813265) q[2];
rz(0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(-0.65778041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0021129) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(-2.9667903) q[0];
rz(2.8305956) q[2];
sx q[2];
rz(-1.9823091) q[2];
sx q[2];
rz(-0.81931537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8448338) q[1];
sx q[1];
rz(-1.3169603) q[1];
sx q[1];
rz(1.0253419) q[1];
rz(-0.25810453) q[3];
sx q[3];
rz(-1.8542284) q[3];
sx q[3];
rz(-0.62844814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(1.0470541) q[2];
rz(-2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.8803546) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(1.0239333) q[0];
rz(-2.3539885) q[1];
sx q[1];
rz(-1.5757631) q[1];
sx q[1];
rz(-0.62686282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1889362) q[0];
sx q[0];
rz(-0.034010012) q[0];
sx q[0];
rz(-2.0859499) q[0];
rz(-1.2704029) q[2];
sx q[2];
rz(-1.2671748) q[2];
sx q[2];
rz(-2.6189569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31506854) q[1];
sx q[1];
rz(-1.8693923) q[1];
sx q[1];
rz(1.8930757) q[1];
x q[2];
rz(1.2792148) q[3];
sx q[3];
rz(-0.40816187) q[3];
sx q[3];
rz(-0.51418958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(-2.5300238) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(1.942379) q[0];
rz(1.9723643) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(2.8170524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31774662) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(-0.86156396) q[0];
rz(-2.3338823) q[2];
sx q[2];
rz(-1.0683904) q[2];
sx q[2];
rz(0.051256996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.868184) q[1];
sx q[1];
rz(-1.5955828) q[1];
sx q[1];
rz(-2.1048057) q[1];
rz(-pi) q[2];
rz(0.59431521) q[3];
sx q[3];
rz(-2.6404877) q[3];
sx q[3];
rz(1.8584115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67363182) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(-1.8661873) q[2];
rz(-2.0868789) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6773029) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(-1.7328847) q[0];
rz(-0.45577058) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.9394402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5613272) q[0];
sx q[0];
rz(-0.94764493) q[0];
sx q[0];
rz(-2.6448963) q[0];
rz(-pi) q[1];
rz(-1.1832778) q[2];
sx q[2];
rz(-2.2829208) q[2];
sx q[2];
rz(2.9930263) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0899883) q[1];
sx q[1];
rz(-3.0172536) q[1];
sx q[1];
rz(0.73595562) q[1];
rz(-pi) q[2];
rz(-0.2644516) q[3];
sx q[3];
rz(-1.3925941) q[3];
sx q[3];
rz(1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1414286) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.6092469) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5672011) q[0];
sx q[0];
rz(-0.53744692) q[0];
sx q[0];
rz(-2.6525932) q[0];
rz(-2.3959827) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(1.9288043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6541877) q[1];
sx q[1];
rz(-1.1206756) q[1];
sx q[1];
rz(0.022189157) q[1];
x q[2];
rz(0.17144449) q[3];
sx q[3];
rz(-2.532634) q[3];
sx q[3];
rz(-2.3122548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9939076) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.9130075) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.6555697) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-2.4771931) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(-1.2845576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212095) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(2.3032805) q[0];
rz(2.1696854) q[2];
sx q[2];
rz(-2.6910057) q[2];
sx q[2];
rz(2.633) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6397275) q[1];
sx q[1];
rz(-2.0149391) q[1];
sx q[1];
rz(0.72413866) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3900208) q[3];
sx q[3];
rz(-1.5454925) q[3];
sx q[3];
rz(-2.7109773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61218843) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-2.5908296) q[2];
rz(-0.70358706) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7228912) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89833242) q[0];
sx q[0];
rz(-1.4295477) q[0];
sx q[0];
rz(-0.54153533) q[0];
rz(2.6780307) q[2];
sx q[2];
rz(-2.2710685) q[2];
sx q[2];
rz(2.314032) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46882792) q[3];
sx q[3];
rz(-0.83823293) q[3];
sx q[3];
rz(1.9459141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71173944) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(1.343887) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(0.96991878) q[2];
sx q[2];
rz(-0.51545943) q[2];
sx q[2];
rz(-3.1209844) q[2];
rz(-2.5034954) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
