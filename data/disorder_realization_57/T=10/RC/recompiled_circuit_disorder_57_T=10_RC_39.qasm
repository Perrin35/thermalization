OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(2.6101987) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019842) q[0];
sx q[0];
rz(-1.3586205) q[0];
sx q[0];
rz(-0.19594812) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4389624) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(-1.4544912) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0002999) q[1];
sx q[1];
rz(-1.3360026) q[1];
sx q[1];
rz(1.512212) q[1];
x q[2];
rz(-0.12227998) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(1.6137326) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(-2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(-0.81940991) q[0];
rz(2.8557414) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.2664638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0985384) q[0];
sx q[0];
rz(-2.1715954) q[0];
sx q[0];
rz(-1.1061125) q[0];
rz(1.9752713) q[2];
sx q[2];
rz(-1.7911439) q[2];
sx q[2];
rz(-1.5420367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5310865) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(0.83341877) q[1];
rz(-2.5655079) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(-1.2036122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1566029) q[2];
sx q[2];
rz(-1.9282324) q[2];
sx q[2];
rz(1.161969) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(0.073908977) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6067628) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8436463) q[0];
sx q[0];
rz(-2.2046411) q[0];
sx q[0];
rz(-1.6266536) q[0];
rz(-pi) q[1];
rz(-0.60909231) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(2.7012205) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68413098) q[1];
sx q[1];
rz(-2.0254454) q[1];
sx q[1];
rz(0.40747868) q[1];
rz(0.89214274) q[3];
sx q[3];
rz(-2.3193079) q[3];
sx q[3];
rz(-0.014978623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(1.5807318) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(2.4981807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075539) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(0.87189829) q[0];
rz(2.7543228) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(-0.29104582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5975808) q[0];
sx q[0];
rz(-2.512092) q[0];
sx q[0];
rz(1.175571) q[0];
rz(1.8918858) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(0.89821399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93451553) q[1];
sx q[1];
rz(-1.0672489) q[1];
sx q[1];
rz(-2.3758923) q[1];
rz(1.325874) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(0.50928515) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0020224) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(0.64506662) q[0];
rz(0.6257261) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-2.5508946) q[0];
sx q[0];
rz(0.22038711) q[0];
rz(-1.9899881) q[2];
sx q[2];
rz(-0.79608166) q[2];
sx q[2];
rz(-2.5615356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1944151) q[1];
sx q[1];
rz(-1.0395323) q[1];
sx q[1];
rz(0.26015014) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.418525) q[3];
sx q[3];
rz(-1.0183183) q[3];
sx q[3];
rz(0.0020364062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(1.7774263) q[2];
rz(1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.0548148) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(2.3278918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700304) q[0];
sx q[0];
rz(-0.22261482) q[0];
sx q[0];
rz(1.2827669) q[0];
x q[1];
rz(-2.5021624) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(-2.8384428) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54436436) q[1];
sx q[1];
rz(-1.9042943) q[1];
sx q[1];
rz(-2.3436105) q[1];
rz(-2.2807062) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-0.83356726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(2.4555092) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(0.66326052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2880741) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(3.1115565) q[0];
x q[1];
rz(-1.3528321) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(-1.8916212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.315553) q[1];
sx q[1];
rz(-1.3693046) q[1];
sx q[1];
rz(1.0793346) q[1];
rz(-1.0955986) q[3];
sx q[3];
rz(-2.0810658) q[3];
sx q[3];
rz(2.1197532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-3.1265124) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(2.6838141) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37516884) q[0];
sx q[0];
rz(-0.94978118) q[0];
sx q[0];
rz(-0.39874886) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56577487) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(2.4772252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.01955186) q[1];
sx q[1];
rz(-1.3380088) q[1];
sx q[1];
rz(-1.9368534) q[1];
x q[2];
rz(0.55191603) q[3];
sx q[3];
rz(-1.8822) q[3];
sx q[3];
rz(-2.6375254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(2.5568533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027771587) q[0];
sx q[0];
rz(-0.50243938) q[0];
sx q[0];
rz(2.9071945) q[0];
x q[1];
rz(1.2051177) q[2];
sx q[2];
rz(-2.0846016) q[2];
sx q[2];
rz(-2.6363381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4084088) q[1];
sx q[1];
rz(-1.6400664) q[1];
sx q[1];
rz(-2.7015711) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0869914) q[3];
sx q[3];
rz(-0.30011794) q[3];
sx q[3];
rz(2.8407612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(-2.8653223) q[2];
rz(-2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(-2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722512) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-2.4861091) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(-1.6190593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(-1.9615016) q[0];
rz(-0.68799512) q[2];
sx q[2];
rz(-1.3032459) q[2];
sx q[2];
rz(-0.5984532) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5787705) q[1];
sx q[1];
rz(-0.65014231) q[1];
sx q[1];
rz(-0.27362089) q[1];
rz(-pi) q[2];
rz(-0.75928648) q[3];
sx q[3];
rz(-0.55428234) q[3];
sx q[3];
rz(2.5169773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9860501) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(2.299451) q[2];
rz(1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(3.1150637) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90606541) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(1.7799829) q[2];
sx q[2];
rz(-1.5968512) q[2];
sx q[2];
rz(-3.130393) q[2];
rz(3.0620861) q[3];
sx q[3];
rz(-2.2600997) q[3];
sx q[3];
rz(-3.1117677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
