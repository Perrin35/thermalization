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
rz(1.4580392) q[0];
sx q[0];
rz(5.3661348) q[0];
sx q[0];
rz(9.6714749) q[0];
rz(3.9858272) q[1];
sx q[1];
rz(1.2935473) q[1];
sx q[1];
rz(10.727439) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58606746) q[0];
sx q[0];
rz(-1.2159776) q[0];
sx q[0];
rz(-2.8321502) q[0];
x q[1];
rz(2.6370722) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(2.1736886) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0820391) q[1];
sx q[1];
rz(-2.5055725) q[1];
sx q[1];
rz(-0.18853699) q[1];
rz(2.2095895) q[3];
sx q[3];
rz(-0.18086704) q[3];
sx q[3];
rz(0.34227926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(-1.1945266) q[2];
rz(-1.901769) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164923) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(2.4745353) q[0];
rz(0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(-1.0557231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0843029) q[0];
sx q[0];
rz(-0.94883666) q[0];
sx q[0];
rz(0.33226407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6237359) q[2];
sx q[2];
rz(-2.756697) q[2];
sx q[2];
rz(3.0786849) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42506524) q[1];
sx q[1];
rz(-0.34913539) q[1];
sx q[1];
rz(1.8340111) q[1];
rz(-pi) q[2];
rz(2.3930644) q[3];
sx q[3];
rz(-2.1338042) q[3];
sx q[3];
rz(2.9070118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51862741) q[2];
sx q[2];
rz(-2.3845606) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(2.8957497) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.2639379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-0.73391947) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(2.6626185) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(-1.3165547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265357) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(-0.75124426) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20256217) q[2];
sx q[2];
rz(-1.1656467) q[2];
sx q[2];
rz(0.7140401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90637246) q[1];
sx q[1];
rz(-2.2487469) q[1];
sx q[1];
rz(1.9560019) q[1];
x q[2];
rz(2.1785424) q[3];
sx q[3];
rz(-0.96745771) q[3];
sx q[3];
rz(-2.778229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58671826) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(-0.47373104) q[2];
rz(-0.50940618) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(1.9853076) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(-1.0152869) q[0];
rz(0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-0.36605787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81824077) q[0];
sx q[0];
rz(-2.2466772) q[0];
sx q[0];
rz(0.9267207) q[0];
rz(-pi) q[1];
rz(-0.54738657) q[2];
sx q[2];
rz(-0.27392188) q[2];
sx q[2];
rz(-0.78597921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3234933) q[1];
sx q[1];
rz(-2.202569) q[1];
sx q[1];
rz(-1.1786103) q[1];
x q[2];
rz(0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4801243) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(0.69501957) q[2];
rz(2.2906637) q[3];
sx q[3];
rz(-1.3635819) q[3];
sx q[3];
rz(1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.618914) q[1];
sx q[1];
rz(-0.42627898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.016092) q[0];
sx q[0];
rz(-2.3707485) q[0];
sx q[0];
rz(-1.9941814) q[0];
rz(-pi) q[1];
rz(2.2441909) q[2];
sx q[2];
rz(-0.36172141) q[2];
sx q[2];
rz(2.7182686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6719088) q[1];
sx q[1];
rz(-1.6571181) q[1];
sx q[1];
rz(1.1344019) q[1];
rz(-1.0067389) q[3];
sx q[3];
rz(-2.5814179) q[3];
sx q[3];
rz(-0.87040802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8729426) q[2];
sx q[2];
rz(-0.23317569) q[2];
sx q[2];
rz(-3.1128913) q[2];
rz(-0.21102333) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(0.15897121) q[0];
rz(-0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(1.5370625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814674) q[0];
sx q[0];
rz(-2.1985538) q[0];
sx q[0];
rz(2.5808236) q[0];
rz(2.3889306) q[2];
sx q[2];
rz(-0.48303451) q[2];
sx q[2];
rz(-0.31980425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83035806) q[1];
sx q[1];
rz(-1.0060961) q[1];
sx q[1];
rz(0.021685251) q[1];
x q[2];
rz(-1.2605242) q[3];
sx q[3];
rz(-1.2168443) q[3];
sx q[3];
rz(1.1502532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0663466) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510821) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(2.2313927) q[0];
rz(-0.74288145) q[1];
sx q[1];
rz(-1.0123092) q[1];
sx q[1];
rz(1.9076617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54523477) q[0];
sx q[0];
rz(-1.6420134) q[0];
sx q[0];
rz(-0.13990732) q[0];
rz(-pi) q[1];
rz(-1.9803488) q[2];
sx q[2];
rz(-1.8170333) q[2];
sx q[2];
rz(-2.4135532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90335315) q[1];
sx q[1];
rz(-0.80447865) q[1];
sx q[1];
rz(2.5969905) q[1];
rz(-pi) q[2];
rz(1.1804092) q[3];
sx q[3];
rz(-1.015706) q[3];
sx q[3];
rz(0.11484717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-2.2255955) q[2];
rz(0.25137526) q[3];
sx q[3];
rz(-2.1706457) q[3];
sx q[3];
rz(2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.4995572) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(-3.1009951) q[0];
rz(-1.9675072) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(2.0814799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745251) q[0];
sx q[0];
rz(-0.39345783) q[0];
sx q[0];
rz(0.25263806) q[0];
rz(-pi) q[1];
rz(1.6868851) q[2];
sx q[2];
rz(-1.5512067) q[2];
sx q[2];
rz(1.7766812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1402982) q[1];
sx q[1];
rz(-2.4457275) q[1];
sx q[1];
rz(0.49927478) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4701263) q[3];
sx q[3];
rz(-1.1219624) q[3];
sx q[3];
rz(-1.1374813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(2.3504284) q[2];
rz(1.5353954) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(-2.9839363) q[0];
rz(-3.0158896) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(0.30473614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19491542) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(-1.2437245) q[0];
rz(-2.5013148) q[2];
sx q[2];
rz(-1.9140179) q[2];
sx q[2];
rz(-2.8514112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.803639) q[1];
sx q[1];
rz(-1.5731166) q[1];
sx q[1];
rz(0.13428899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7823506) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(-0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2711082) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(1.6893207) q[2];
rz(2.3952386) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8993503) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(0.53264701) q[0];
rz(-0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(-2.9740082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78412752) q[0];
sx q[0];
rz(-1.7622158) q[0];
sx q[0];
rz(-3.0166059) q[0];
x q[1];
rz(2.8791076) q[2];
sx q[2];
rz(-1.0392611) q[2];
sx q[2];
rz(2.558311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1780562) q[1];
sx q[1];
rz(-1.85317) q[1];
sx q[1];
rz(-2.8368188) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3066809) q[3];
sx q[3];
rz(-1.2560227) q[3];
sx q[3];
rz(-1.9251458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.31183895) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(2.6746542) q[2];
rz(-2.0385108) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(-1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5760096) q[0];
sx q[0];
rz(-0.82721114) q[0];
sx q[0];
rz(0.4785434) q[0];
rz(-0.78454984) q[1];
sx q[1];
rz(-0.34930925) q[1];
sx q[1];
rz(0.92225155) q[1];
rz(-0.68741531) q[2];
sx q[2];
rz(-0.981642) q[2];
sx q[2];
rz(-2.8264075) q[2];
rz(0.26950027) q[3];
sx q[3];
rz(-1.2242347) q[3];
sx q[3];
rz(0.4927529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
