OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(2.9076599) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(4.5555727) q[1];
sx q[1];
rz(11.04784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9568125) q[0];
sx q[0];
rz(-2.0818458) q[0];
sx q[0];
rz(1.2132116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1351003) q[2];
sx q[2];
rz(-2.6681136) q[2];
sx q[2];
rz(1.651929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17659345) q[1];
sx q[1];
rz(-1.0142583) q[1];
sx q[1];
rz(2.2189369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99125864) q[3];
sx q[3];
rz(-2.391444) q[3];
sx q[3];
rz(-1.0446435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(-0.8055299) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.7854179) q[3];
sx q[3];
rz(-0.75683561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(0.43847325) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(-3.0335887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9697125) q[0];
sx q[0];
rz(-1.6816655) q[0];
sx q[0];
rz(3.1165439) q[0];
x q[1];
rz(-2.503452) q[2];
sx q[2];
rz(-2.2109988) q[2];
sx q[2];
rz(-1.3253044) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9671772) q[1];
sx q[1];
rz(-1.4870054) q[1];
sx q[1];
rz(1.3166974) q[1];
rz(1.2233954) q[3];
sx q[3];
rz(-0.51651556) q[3];
sx q[3];
rz(-1.4054738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3652304) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(2.727437) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(2.3057002) q[0];
rz(1.6235141) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(1.2380884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94592124) q[0];
sx q[0];
rz(-0.64650457) q[0];
sx q[0];
rz(2.2158428) q[0];
x q[1];
rz(0.34332163) q[2];
sx q[2];
rz(-1.3385217) q[2];
sx q[2];
rz(-0.49988036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39565428) q[1];
sx q[1];
rz(-1.6881583) q[1];
sx q[1];
rz(3.0836041) q[1];
rz(-pi) q[2];
rz(-0.79583056) q[3];
sx q[3];
rz(-2.748877) q[3];
sx q[3];
rz(-2.8098729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28430024) q[2];
sx q[2];
rz(-0.14480536) q[2];
sx q[2];
rz(2.4227552) q[2];
rz(-2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6179287) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(-2.0148328) q[0];
rz(-1.2976546) q[1];
sx q[1];
rz(-1.6512066) q[1];
sx q[1];
rz(-2.0984971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27985036) q[0];
sx q[0];
rz(-1.540483) q[0];
sx q[0];
rz(1.8413196) q[0];
rz(0.2506855) q[2];
sx q[2];
rz(-2.0935241) q[2];
sx q[2];
rz(-1.2831068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35395393) q[1];
sx q[1];
rz(-1.2661722) q[1];
sx q[1];
rz(2.980176) q[1];
x q[2];
rz(0.65535069) q[3];
sx q[3];
rz(-1.5028302) q[3];
sx q[3];
rz(-0.36253906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8635233) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-1.0370022) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-1.0629531) q[3];
sx q[3];
rz(0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0015513) q[0];
sx q[0];
rz(-2.934444) q[0];
sx q[0];
rz(-0.087604372) q[0];
rz(-2.7032848) q[1];
sx q[1];
rz(-1.0821082) q[1];
sx q[1];
rz(-1.3137438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739819) q[0];
sx q[0];
rz(-1.0896519) q[0];
sx q[0];
rz(0.26480459) q[0];
rz(-pi) q[1];
rz(-1.472166) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(-1.94869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96927724) q[1];
sx q[1];
rz(-1.7427497) q[1];
sx q[1];
rz(0.090222619) q[1];
rz(-0.99788061) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(-2.5469766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(-0.56337774) q[2];
rz(1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(2.5105072) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1221984) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(0.011938183) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48680533) q[0];
sx q[0];
rz(-1.1410895) q[0];
sx q[0];
rz(-0.91581099) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8831364) q[2];
sx q[2];
rz(-1.788013) q[2];
sx q[2];
rz(-2.1803792) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4820741) q[1];
sx q[1];
rz(-2.4122752) q[1];
sx q[1];
rz(-0.68406711) q[1];
rz(-2.7390476) q[3];
sx q[3];
rz(-1.8410826) q[3];
sx q[3];
rz(0.2337993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7616854) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948792) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(0.43103257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1566136) q[0];
sx q[0];
rz(-2.7811858) q[0];
sx q[0];
rz(-0.63215881) q[0];
x q[1];
rz(0.22054976) q[2];
sx q[2];
rz(-0.23633453) q[2];
sx q[2];
rz(-2.5823809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7546852) q[1];
sx q[1];
rz(-2.0758938) q[1];
sx q[1];
rz(-2.7688857) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0484974) q[3];
sx q[3];
rz(-2.5601031) q[3];
sx q[3];
rz(3.1149816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3057574) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(0.26710278) q[2];
rz(3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-3.035868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(-2.5015976) q[0];
rz(-0.85482875) q[1];
sx q[1];
rz(-1.6565485) q[1];
sx q[1];
rz(1.7291501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4025637) q[0];
sx q[0];
rz(-1.1629346) q[0];
sx q[0];
rz(0.5824851) q[0];
rz(-pi) q[1];
rz(-1.475263) q[2];
sx q[2];
rz(-2.6474806) q[2];
sx q[2];
rz(-1.4612559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58929491) q[1];
sx q[1];
rz(-1.4865685) q[1];
sx q[1];
rz(0.34119795) q[1];
rz(-pi) q[2];
rz(-2.2858587) q[3];
sx q[3];
rz(-0.54045709) q[3];
sx q[3];
rz(3.1202701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60014805) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(-2.6680434) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(0.1951018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.37050978) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(-1.3516082) q[0];
rz(0.31556684) q[1];
sx q[1];
rz(-1.6866997) q[1];
sx q[1];
rz(-1.9884761) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0457927) q[0];
sx q[0];
rz(-1.519916) q[0];
sx q[0];
rz(-1.6232383) q[0];
rz(-pi) q[1];
rz(-1.6193689) q[2];
sx q[2];
rz(-2.1056386) q[2];
sx q[2];
rz(0.64316197) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0259195) q[1];
sx q[1];
rz(-0.5957091) q[1];
sx q[1];
rz(-2.9133948) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38001506) q[3];
sx q[3];
rz(-1.9560062) q[3];
sx q[3];
rz(2.4917701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(-2.7678164) q[2];
rz(-1.2260381) q[3];
sx q[3];
rz(-0.91807476) q[3];
sx q[3];
rz(1.3396243) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.8208338) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(2.6722867) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8784638) q[0];
sx q[0];
rz(-1.6471905) q[0];
sx q[0];
rz(-1.6883649) q[0];
rz(-pi) q[1];
rz(-1.2861757) q[2];
sx q[2];
rz(-0.60065833) q[2];
sx q[2];
rz(0.17135581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9633858) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(0.15060138) q[1];
rz(0.68113459) q[3];
sx q[3];
rz(-0.38060846) q[3];
sx q[3];
rz(-1.8073624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70539537) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(1.2104642) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(-2.1657522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198915) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(-0.89314356) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(-2.2371815) q[2];
sx q[2];
rz(-2.7795962) q[2];
sx q[2];
rz(-2.4301651) q[2];
rz(-1.0004956) q[3];
sx q[3];
rz(-1.199493) q[3];
sx q[3];
rz(-1.5293157) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
