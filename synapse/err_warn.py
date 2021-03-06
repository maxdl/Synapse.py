import sys


class ProfileError(Exception):
    def __init__(self, profile, msg):
        self.profile = profile
        self.msg = msg + "."


def profile_warning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True


def profile_message(msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)